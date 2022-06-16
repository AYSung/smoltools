"""Functions for selecting residues and atoms from PDB structure."""

from itertools import chain

from Bio.PDB.Atom import Atom
from Bio.PDB.Chain import Chain
from Bio.PDB.Entity import Entity
from Bio.PDB.Residue import Residue
from Bio.PDB import Selection
from Bio.PDB.Structure import Structure


def get_chain(structure: Structure, model: int = 0, chain: str = 'A') -> Chain:
    """
    Returns a chain from a PDB structure object. Default is to return chain A of the
    first model.

    Parameters:
    -----------
    structure (Structure): PDB structure object.
    model (int): Model number, default is 0.
    chain (str): Chain identifier, default is A.

    Returns:
    --------
    Chain: PDB chain object.
    """
    return structure[model][chain]


def get_residues(entity: Entity, residue_filter: set[str] = None) -> list[Residue]:
    """
    Produces a list of all residues in a PDB entity. Can provide a set of specific
        residues to keep.

    Parameters:
    -----------
    entity (Entity): PDB entity object.
    residue_filter (set[str]): Optional, a set (or other list-like) of three letter
        amino codes for the residues to keep. Default is to return all residues.

    Returns:
    --------
    list[Residue]: List of PDB residue objects in the given entity that meet the
        residue filter.
    """
    residues = Selection.unfold_entities(entity, 'R')
    if residue_filter is None:
        return residues
    else:
        return list(filter(lambda x: x.get_resname() in residue_filter, residues))


def get_alpha_carbons(residues: list[Residue]) -> list[Atom]:
    """
    Returns a list of alpha carbons for a given list of residues.

    Parameters:
    -----------
    residues (list[Residue]): list of PDB residue objects.

    Returns:
    --------
    list[Atom]: list of alpha carbons as PDB atom objects.
    """
    return [residue['CA'] for residue in residues]


def get_carbons(
    residues: list[Residue], atom_select: dict[str : list[str]]
) -> list[Atom]:
    """
    Returns a list of atoms from a list of residues that meet the atom selection
    criteria. Requires a dictionary of the names of the atoms to retrieve for each
    amino acid.
    """

    def _get_atoms(residue: Residue, atom_names: list[str]):
        return [residue[atom_name] for atom_name in atom_names]

    # def _filter_atom_name(residue: Residue) -> list[Atom]:
    #     """ """
    #     return list(
    #         filter(
    #             lambda atom: atom.get_name() in atom_select[residue.get_resname()],
    #             residue.get_atoms(),
    #         )
    #     )

    atoms = [
        _get_atoms(residue, atom_select[residue.get_resname()]) for residue in residues
    ]
    return list(chain(*atoms))
