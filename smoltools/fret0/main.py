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
