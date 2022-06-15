from Bio.PDB.Structure import Structure
import pandas as pd

from fret0.distance import pairwise_distance_between_conformations
from fret0.pdb_parser import read_atom_coords, get_sasa
from fret0.residue_filter import (
    combine_surface_residues,
    filter_pairwise_distances_by_sasa,
)


def pairwise_distances(
    structure_a: Structure,
    chain_a: str,
    structure_b: Structure,
    chain_b: str,
    sasa_cutoff: float = None,
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
    coords_a = read_atom_coords(structure_a, chain_a)
    coords_b = read_atom_coords(structure_b, chain_b)
    pairwise_distance_df = pairwise_distance_between_conformations(coords_a, coords_b)

    if sasa_cutoff is not None:
        sasa_a = get_sasa(structure_a, chain_a)
        sasa_b = get_sasa(structure_a, chain_b)
        surface_residues = combine_surface_residues(sasa_a, sasa_b, cutoff=sasa_cutoff)
        pairwise_distance_df = filter_pairwise_distances_by_sasa(
            pairwise_distance_df, surface_residues
        )

    return pairwise_distance_df
