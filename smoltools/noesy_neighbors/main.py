from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
import pandas as pd

from smoltools.pdbtools import path_to_chain, coordinate_table
import smoltools.pdbtools.select as select


LABELED_CARBONS = {
    'ILV': {
        'ILE': ['CD', 'CD1', 'CG2'],
        'LEU': ['CD1', 'CD2'],
        'VAL': ['CG1', 'CG2'],
    },
    'ILVA': {
        'ILE': ['CD', 'CD1'],
        'LEU': ['CD1', 'CD2'],
        'VAL': ['CG1', 'CG2'],
        'ALA': ['CB'],
    },
    'ILVMAT': {
        'ILE': ['CD', 'CD1', 'CG2'],
        'LEU': ['CD1', 'CD2'],
        'VAL': ['CG1', 'CG2'],
        'MET': ['CE'],
        'ALA': ['CB'],
        'THR': ['CG2'],
    },
}

LABELING_SCHEMES = list(LABELED_CARBONS.keys())


def get_labeled_carbons(
    residues: list[Residue], labeled_atoms: dict[str, list[str]]
) -> list[Atom]:
    """Retrieve labelled carbons from branched-chain amino acids (VAL, LEU, ILE)
    from a list of residues.

    Parameters:
    -----------
    residues (list[Residue]): List of PDB Residue objects.
    labeled_atoms (dict): Dictionary mapping three letter residue ID (e.g. 'ILE')
        to list of atoms to select (e.g. ['CD', 'CG2'])

    Returns:
    list[Atom]: List of PDB Atom objects.
    """

    return select.get_carbons(residues, labeled_atoms)


def coordinates_from_chain(
    chain: Chain, labeled_atoms: dict[str, list[str]]
) -> pd.DataFrame:
    """Calculate pairwise distances of terminal carbons of branched-chain amino acids
    in the given Chain object. Use if a chain object is already loaded.

    Parameters:
    -----------
    chain (Chain): PDB Chain object.
    labeled_atoms (dict): Dictionary mapping three letter residue ID (e.g. 'ILE')
        to list of atoms to select (e.g. ['CD', 'CG2'])

    Returns:
    --------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.
    """
    residue_filter = set(labeled_atoms.keys())
    residues = select.get_residues(chain, residue_filter=residue_filter)
    atoms = get_labeled_carbons(residues, labeled_atoms)
    return (
        coordinate_table(atoms)
        .assign(
            id=lambda x: x.residue_name + x.residue_number.astype(str) + '-' + x.atom_id
        )
        .set_index('id')
        .loc[:, ['x', 'y', 'z']]
    )


def coordinates_from_path(
    path: str,
    labeled_atoms: dict[str, list[str]],
    model: int = 0,
    chain: str = 'A',
) -> pd.DataFrame:
    """Calculate pairwise distances of terminal carbons of branched-chain amino acids
    in the specified chain from a PDB file. Use if starting directly from PDB file.

    Parameters:
    -----------
    path (str): Path to PDB file.
    labeled_atoms (dict): Dictionary mapping three letter residue ID (e.g. 'ILE')
        to list of atoms to select (e.g. ['CD', 'CG2'])
    model (int): Model number of desired chain (default = 0)
    chain (str): Chain ID of desired chain (default = 'A')

    Returns:
    --------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.
    """
    chain = path_to_chain(path, model=model, chain=chain)
    return coordinates_from_chain(chain, labeled_atoms)


def coordinates_from_path_presets(
    path: str,
    mode: str = 'ILV',
    model: int = 0,
    chain: str = 'A',
) -> pd.DataFrame:
    """Calculate pairwise distances of terminal carbons of branched-chain amino acids
    in the specified chain from a PDB file. Use if starting directly from PDB file.

    Parameters:
    -----------
    path (str): Path to PDB file.
    mode (str): Predefined labeled atom selections (choices are 'ILV', 'ILVA', and 'ILVMAT')
    model (int): Model number of desired chain (default = 0)
    chain (str): Chain ID of desired chain (default = 'A')

    Returns:
    --------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.
    """
    labeled_atoms = LABELED_CARBONS[mode]
    chain = path_to_chain(path, model=model, chain=chain)
    return coordinates_from_chain(chain, labeled_atoms)
