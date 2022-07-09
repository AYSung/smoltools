"""Functions for selecting residues and atoms from PDB structure."""

from itertools import chain

from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure

from smoltools.pdbtools.exceptions import ChainNotFound


def get_chain(structure: Structure, model: int, chain: str) -> Chain:
    """Returns a chain from a PDB structure object.

    Parameters:
    -----------
    structure (Structure): PDB structure object.
    model (int): Model number.
    chain (str): Chain identifier.

    Returns:
    --------
    Chain: PDB chain object.
    """
    try:
        return structure[model][chain]
    except KeyError as e:
        raise ChainNotFound(structure.get_id(), model, chain) from e


def get_residues(chain: Chain, residue_filter: set[str] = None) -> list[Residue]:
    """Produces a list of all residues in a PDB chain. Can provide a set of specific
    residues to keep.

    Parameters:
    -----------
    chain (Chain): PDB chain object.
    residue_filter (set[str]): Optional, a set (or other list-like) of three letter
        amino codes for the residues to keep. Default is to return all residues.

    Returns:
    --------
    list[Residue]: List of PDB residue objects in the given entity that meet the
        residue filter.
    """
    residues = [residue for residue in chain.get_residues()]
    if residue_filter is None:
        return [residue for residue in residues if residue.get_id()[0] == ' ']
    else:
        return [
            residue for residue in residues if residue.get_resname() in residue_filter
        ]
    # TODO: error handling for empty residue list?


def get_alpha_carbons(residues: list[Residue]) -> list[Atom]:
    """Returns a list of alpha carbons for a given list of residues.

    Parameters:
    -----------
    residues (list[Residue]): list of PDB residue objects.

    Returns:
    --------
    list[Atom]: list of alpha carbons as PDB atom objects.
    """

    def _get_alpha_carbon(residue: Residue) -> Atom:
        alpha_carbons = [atom for atom in residue.get_atoms() if atom.get_id() == 'CA']
        if not alpha_carbons:
            print(f'{residue} has no alpha carbon')

        return alpha_carbons

    alpha_carbons = [_get_alpha_carbon(residue) for residue in residues]
    return list(chain(*alpha_carbons))
    # TODO: error handling for empty residue list?


def get_carbons(
    residues: list[Residue], atom_select: dict[str : list[str]]
) -> list[Atom]:
    """Returns a list of atoms from a list of residues that meet the atom selection
    criteria. Requires a dictionary of the names of the atoms to retrieve for each
    amino acid.
    """

    def _get_atoms(residue: Residue, atom_names: list[str]):
        return [atom for atom in residue.get_atoms() if atom.get_name() in atom_names]

    atoms = [
        _get_atoms(residue, atom_select[residue.get_resname()]) for residue in residues
    ]
    return list(chain(*atoms))
    # TODO: error handling for empty residue list?


def filter_by_b_factor(atoms: list[Atom], cutoff) -> list[Atom]:
    """Returns a list of atoms with a b factor that meets the provided cutoff."""
    return [atom for atom in atoms if atom.get_bfactor() > cutoff]
