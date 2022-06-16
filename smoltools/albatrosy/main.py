from pathlib import Path

from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue

# import numpy as np
import pandas as pd

import smoltools.calculate.distance as distance
import smoltools.pdbtools.load as load
import smoltools.pdbtools.select as select


def get_labelled_carbons(residues: list[Residue]) -> list[Atom]:
    LABELLED_CARBONS = {
        'VAL': ['CG1', 'CG2'],
        'LEU': ['CD1', 'CD2'],
        'ILE': ['CD1'],
    }

    return select.get_carbons(residues, LABELLED_CARBONS)


def coordinate_table(atoms: list[Atom]) -> pd.DataFrame:
    def _get_atom_info(atom: Atom) -> tuple:
        residue_number = atom.get_parent().get_id()[1]
        return f'{residue_number}-{atom.get_name()}', *atom.get_coord()

    atom_info = [_get_atom_info(atom) for atom in atoms]

    return pd.DataFrame(atom_info, columns=['atom_id', 'x', 'y', 'z']).set_index(
        'atom_id'
    )


def structure_to_distances(
    path: str | Path, model: int = 0, chain: str = 'A'
) -> pd.DataFrame:
    structure = load.read_pdb_from_path(path)
    chain = select.get_chain(structure, model=model, chain=chain)
    residues = select.get_residues(chain, residue_filter={'VAL', 'LEU', 'ILE'})
    labelled_atoms = get_labelled_carbons(residues)
    coords = coordinate_table(labelled_atoms)
    return distance.calculate_pairwise_distances(coords)
