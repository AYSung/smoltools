"""Functions for calculating FRET efficiencies."""

import numpy as np
import pandas as pd


def _calculate_e_fret(distance: pd.Series, r0: float) -> pd.Series:
    """Calculate FRET efficiency based on inter-residue distances for a given r0."""
    return 1 / (1 + (distance / r0) ** 6)


def _calculate_delta_e_fret(e_fret_a: pd.Series, e_fret_b: pd.Series) -> pd.Series:
    """Calculate the magnitude of the difference in FRET efficiencies."""
    return e_fret_a - e_fret_b


def e_fret_between_conformations(df: pd.DataFrame, r0: float) -> pd.DataFrame:
    """Calculate FRET efficiencies from a pairwise distance DataFrame.

    Parameters:
    -----------
    df (DataFrame): DataFrame with pairwise distances for conformation A and
        conformation B.
    r0 (float): R0 values used for calculating FRET efficiency.

    Returns:
    --------
    DataFrame: DataFrame with FRET efficiency calculate for each residue
        pair, as well as the change in FRET efficiency between conformations.
    """
    return df[['id_1', 'id_2']].assign(
        E_fret_a=_calculate_e_fret(df.distance_a, r0),
        E_fret_b=_calculate_e_fret(df.distance_b, r0),
        delta_E_fret=lambda x: _calculate_delta_e_fret(x.E_fret_a, x.E_fret_b),
    )


def generate_r0_curve(distance_a: float, distance_b: float) -> pd.DataFrame:
    """Generate data for FRET efficiency as a function of R0 for two distances."""
    r0_range = list(range(20, 81))
    e_fret_a = [_calculate_e_fret(distance_a, r0) for r0 in r0_range]
    e_fret_b = [_calculate_e_fret(distance_b, r0) for r0 in r0_range]
    e_fret = pd.DataFrame(
        {
            'r0': r0_range,
            'A': e_fret_a,
            'B': e_fret_b,
        }
    )

    e_fret_by_distance = e_fret.melt(
        id_vars='r0', var_name='distance', value_name='e_fret'
    )
    e_fret_delta = e_fret.assign(delta=lambda x: np.abs(x.A - x.B)).drop(
        columns=['A', 'B']
    )

    return e_fret_by_distance, e_fret_delta
