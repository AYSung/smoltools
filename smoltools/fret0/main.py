from Bio.PDB.Chain import Chain
import pandas as pd

import smoltools.calculate.distance as distance
from smoltools.pdbtools import path_to_chain, coordinate_table
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
    coords = (
        coordinate_table(alpha_carbons)
        .assign(id=lambda x: x.residue_name + x.residue_number.astype(str))
        .set_index('id')
        .loc[:, ['x', 'y', 'z']]
    )
    return distance.pairwise_distances(coords)


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
