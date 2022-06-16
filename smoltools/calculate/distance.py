"""Functions for calculating atomic distances."""

import numpy as np
import pandas as pd

import scipy.spatial.distance as ssd


def _pairwise_distance(df: pd.DataFrame) -> np.ndarray:
    return ssd.cdist(df, df, 'euclidean')


def _tidy_pairwise_distances(df: pd.DataFrame) -> pd.DataFrame:
    """Take a square dataframe of pairwise distances and convert it to tidy format."""
    return df.melt(value_name='distance', ignore_index=False).reset_index()


def calculate_pairwise_distances(df: pd.DataFrame) -> pd.DataFrame:
    """Given a dataframe with 3D coordinates of each residue, calculate the pairwise
    distance between each residue.
    """
    return (
        pd.DataFrame(
            _pairwise_distance(df),
            index=df.index,
            columns=df.index,
        )
        .rename_axis(index='atom_id_1', columns='atom_id_2')
        .pipe(_tidy_pairwise_distances)
    )


def _merge_pairwise_distances(df_a: pd.DataFrame, df_b: pd.DataFrame) -> pd.DataFrame:
    """Merge two DataFrames of pairwise distances (intersection of residues pairs in each
    dataset)
    """
    return pd.merge(
        df_a,
        df_b,
        on=['atom_id_1', 'atom_id_2'],
        suffixes=['_a', '_b'],
        validate='1:1',
    )


def _calculate_delta_distance(
    distance_a: pd.Series, distance_b: pd.Series
) -> pd.Series:
    """Absolute value of the difference between pairwise differences between the two
    conformations.
    """
    return np.abs(distance_a - distance_b)


def pairwise_distance_between_conformations(
    coord_a: pd.DataFrame, coord_b: pd.DataFrame
) -> pd.DataFrame:
    """
    Given two DataFrames with the 3D coordinates of residues in two conformations,
    return the pairwise distances between each residue in each conformation, as well
    as the difference in pairwise distances between the two conformations.

    Args:
        coord_a (DataFrame): DataFrame with x, y, z coordinates of residues in
            conformation A.
        coord_b (DataFrame): DataFrame with x, y, z coordinates of residues in
            conformation B.

    Returns:
        DataFrame: pandas DataFrame of the pairwise distance between residues for each
            of the two conformations, and the difference in the pairwise distances
            between the conformations. distances reported in angstroms.
    """
    distance_a = calculate_pairwise_distances(coord_a)
    distance_b = calculate_pairwise_distances(coord_b)

    return _merge_pairwise_distances(distance_a, distance_b,).assign(
        delta_distance=lambda x: _calculate_delta_distance(x.distance_a, x.distance_b)
    )
