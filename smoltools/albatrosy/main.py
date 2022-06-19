from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
import pandas as pd

import smoltools.calculate.distance as distance
import smoltools.pdbtools.load as load
import smoltools.pdbtools.select as select


def get_labelled_carbons(residues: list[Residue]) -> list[Atom]:
    """Retrieve labelled carbons from branched-chain amino acids (VAL, LEU, ILE)
    from a list of residues.

    Parameters:
    -----------
    residues (list[Residue]): List of PDB Residue objects.

    Returns:
    list[Atom]: List of PDB Atom objects.
    """
    LABELLED_CARBONS = {
        'VAL': ['CG1', 'CG2'],
        'LEU': ['CD1', 'CD2'],
        'ILE': ['CD1'],
    }

    return select.get_carbons(residues, LABELLED_CARBONS)


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
        return f'{residue_number}-{atom.get_name()}', *atom.get_coord()

    atom_info = [_get_atom_info(atom) for atom in atoms]

    return pd.DataFrame(atom_info, columns=['atom_id', 'x', 'y', 'z']).set_index(
        'atom_id'
    )


def path_to_chain(path: str, model: int = 0, chain: str = 'A') -> Chain:
    """Extract the specified chain from a PDB file.

    Parameters:
    -----------
    path (str): Path to PDB file.
    model (int): Model number of desired chain (default = 0)
    chain (str): Chain ID of desired chain (default = 'A')

    Returns:
    --------
    Chain: PDB Chain object.
    """
    structure = load.read_pdb_from_path(path)
    return select.get_chain(structure, model=model, chain=chain)


def chain_to_distances(chain: Chain) -> pd.DataFrame:
    """Extract 3D coordinates of terminal carbons of branched-chain amino acids in the
    given Chain object. Use if a chain object is already loaded.

    Parameters:
    -----------
    chain (Chain): PDB Chain object.

    Returns:
    --------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.
    """
    residues = select.get_residues(chain, residue_filter={'VAL', 'LEU', 'ILE'})
    labelled_atoms = get_labelled_carbons(residues)
    coords = coordinate_table(labelled_atoms)
    return distance.calculate_pairwise_distances(coords)


def path_to_distances(path: str, model: int = 0, chain: str = 'A') -> pd.DataFrame:
    """Extract 3D coordinates of terminal carbons of branched-chain amino acids in the
    specified chain from a PDB file. Use if starting directly from PDB file.

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
    return chain_to_distances(chain)
