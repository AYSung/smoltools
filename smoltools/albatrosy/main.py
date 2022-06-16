from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

# import numpy as np
import pandas as pd

from smoltools.pdb.select import get_carbons


def get_labelled_carbons(residues: list[Residue]) -> list[Atom]:
    LABELLED_CARBONS = {
        'VAL': ['CG1', 'CG2'],
        'LEU': ['CD1', 'CD2'],
        'ILE': ['CD1'],
    }

    return get_carbons(residues, LABELLED_CARBONS)


def coordinate_table(atoms: list[Atom]) -> pd.DataFrame:
    def _get_atom_info(atom: Atom) -> tuple:
        residue_number = atom.get_parent().get_id()[1]
        return f'{residue_number}-{atom.get_name()}', *atom.get_coord()

    atom_info = [_get_atom_info(atom) for atom in atoms]

    return pd.DataFrame(atom_info, columns=['atom_id', 'x', 'y', 'z']).set_index(
        'atom_id'
    )
