"""Convert list of residues to a table of atomic coordinates"""
from Bio.PDB.Atom import Atom
import pandas as pd


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
        parent_residue = atom.get_parent()
        residue_number = parent_residue.get_id()[1]
        residue_name = parent_residue.get_resname()
        atom_id = atom.get_id()
        coordinates = atom.get_coord()

        return residue_name, residue_number, atom_id, *coordinates

    atom_info = [_get_atom_info(atom) for atom in atoms]
    info_columns = ['residue_name', 'residue_number', 'atom_id']

    return pd.DataFrame(
        atom_info,
        columns=[*info_columns, 'x', 'y', 'z'],
    )
