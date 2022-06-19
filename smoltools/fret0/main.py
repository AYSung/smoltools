from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
import pandas as pd

import smoltools.calculate.distance as distance
from smoltools.pdbtools import path_to_chain
import smoltools.pdbtools.select as select


def chain_to_distances(chain: Chain, sasa_cutoff: float = None) -> pd.DataFrame:
    """Calculate pairwise distances of alpha carbons in the given Chain object.
    Use if a chain object is already loaded.

    Parameters:
    -----------
    chain (Chain): PDB Chain object.

    Returns:
    --------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.
    """
    residues = select.get_residues(chain)
    alpha_carbons = select.get_alpha_carbons(residues)
    if sasa_cutoff is not None:
        alpha_carbons = select.filter_by_b_factor(alpha_carbons, cutoff=sasa_cutoff)
    coords = coordinate_table(alpha_carbons)
    return distance.calculate_pairwise_distances(coords)


def path_to_distances(
    path: str, model: int = 0, chain: str = 'A', sasa_cutoff: float = None
) -> pd.DataFrame:
    """Calculate pairwise distances of alpha carbons in the given Chain object.
    Use if starting directly from PDB file.

    Parameters:
    -----------
    path (str): Path to PDB file.
    model (int): Model number of desired chain (default = 0)
    chain (str): Chain ID of desired chain (default = 'A')

    Returns:
    --------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.
    """
    chain = path_to_chain(path, model=model, chain=chain)
    return chain_to_distances(chain, sasa_cutoff=sasa_cutoff)


def coordinate_table(atoms: list[Atom]) -> pd.DataFrame:
    """Extract 3D coordinates from list of atoms into DataFrame.

    Parameters:
    -----------
    atoms (list[Atom]): List of PDB Atom.

    Returns:
    --------
    DataFrame: Dataframe with the atom ID (residue number, carbon ID) as the index
        and the x, y, z coordinate of each atom as the columns.
    """

    def _get_atom_info(atom: Atom) -> tuple:
        residue_number = atom.get_parent().get_id()[1]
        return residue_number, *atom.get_coord()

    atom_info = [_get_atom_info(atom) for atom in atoms]

    return pd.DataFrame(atom_info, columns=['atom_id', 'x', 'y', 'z']).set_index(
        'atom_id'
    )
