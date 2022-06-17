from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
import pandas as pd

import smoltools.calculate.distance as distance
import smoltools.pdbtools.load as load
import smoltools.pdbtools.select as select


def structure_to_chain(path: str, model: int = 0, chain: str = 'A') -> Chain:
    structure = load.read_pdb_from_path(path)
    return select.get_chain(structure, model=model, chain=chain)


def chain_to_distances(chain: Chain, sasa_cutoff: float = None) -> pd.DataFrame:
    residues = select.get_residues(chain)
    alpha_carbons = select.get_alpha_carbons(residues)
    if sasa_cutoff is not None:
        alpha_carbons = select.filter_by_b_factor(alpha_carbons, cutoff=sasa_cutoff)
    coords = coordinate_table(alpha_carbons)
    return distance.calculate_pairwise_distances(coords)


def structure_to_distances(
    path: str, model: int = 0, chain: str = 'A', sasa_cutoff: float = None
) -> pd.DataFrame:
    chain = structure_to_chain(path, model=model, chain=chain)
    return chain_to_distances(chain, sasa_cutoff=sasa_cutoff)


def coordinate_table(atoms: list[Atom]) -> pd.DataFrame:
    def _get_atom_info(atom: Atom) -> tuple:
        residue_number = atom.get_parent().get_id()[1]
        return residue_number, *atom.get_coord()

    atom_info = [_get_atom_info(atom) for atom in atoms]

    return pd.DataFrame(atom_info, columns=['atom_id', 'x', 'y', 'z']).set_index(
        'atom_id'
    )


def pairwise_distances(
    distances_a: pd.DataFrame,
    distances_b: pd.DataFrame,
) -> pd.DataFrame:
    """
    Reads in the alpha carbon coordinates of two pdb structures representing two
    conformations of the same protein and returns a dataframe of the difference
    in pairwise distances between residues between the two conformations.

    Args:
        structure_a (Structure): pdb Structure object for conformation A
        chain_a (str): chain id for protein of interest in the pdb file for conformation A
        structure_b (Structure): pdb Structure object for conformation B
        chain_b (str): chain id for protein of interest in the pdb file for conformation B
        sasa_cutoff (float): optional value for filtering pairwise distances that only
            contain comparisons between surface exposed residues; requires that the
            relative solvent accessible surface area of each residue be loaded into the
            b factor column of both pdb files.

    Returns:
        DataFrame: pandas DataFrame with the pairwise distances between residues in each
            conformation as well as the difference in pairwise distances between
            conformations.
    """

    return distance._merge_pairwise_distances(distances_a, distances_b).assign(
        delta_distance=lambda x: distance._calculate_delta_distance(
            x.distance_a, x.distance_b
        )
    )
