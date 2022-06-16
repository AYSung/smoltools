import altair as alt
import numpy as np
import pandas as pd

from smoltools.calculate.distance import _merge_pairwise_distances
from smoltools.albatrosy.utils import splice_conformation_tables


def distance_map(df: pd.DataFrame) -> alt.Chart:
    return (
        alt.Chart(df)
        .mark_rect()
        .encode(
            x=alt.X('atom_id_1', title='Atom ID', sort=None),
            y=alt.Y('atom_id_2', title='Atom ID', sort=None),
            color='distance',
            tooltip=[
                alt.Tooltip('atom_id_1', title='Atom #1'),
                alt.Tooltip('atom_id_2', title='Atom #2'),
                alt.Tooltip('distance', title='Distance (\u212B)', format='.1f'),
            ],
        )
        .properties(
            width=600,
            height=600,
        )
    )


def _add_noe_bins(df: pd.DataFrame) -> pd.DataFrame:
    return df.assign(
        noe_strength=lambda x: pd.cut(
            x.distance,
            bins=[0, 2.5, 3.5, 5, 10, np.inf],
            include_lowest=True,
            labels=['strong', 'medium', 'weak', 'very weak', 'none'],
        )
    )


def noe_map(df: pd.DataFrame) -> alt.Chart:
    return (
        alt.Chart(df.pipe(_add_noe_bins))
        .mark_rect()
        .encode(
            x=alt.X('atom_id_1', title='atom ID', sort=None),
            y=alt.Y('atom_id_2', title='atom ID', sort=None),
            color=alt.Color(
                'noe_strength',
                title='NOE',
                scale=alt.Scale(
                    domain=['strong', 'medium', 'weak', 'very weak', 'none'],
                    scheme='blues',
                    reverse=True,
                ),
            ),
            tooltip=[
                alt.Tooltip('atom_id_1', title='Atom #1'),
                alt.Tooltip('atom_id_2', title='Atom #2'),
                alt.Tooltip('distance', title='Distance (\u212B)', format='.1f'),
                alt.Tooltip('noe_strength', title='NOE'),
            ],
        )
        .properties(
            width=600,
            height=600,
        )
    )


def spliced_noe_map(df_a: pd.DataFrame, df_b: pd.DataFrame) -> alt.Chart:
    spliced_df = splice_conformation_tables(df_a, df_b)
    return noe_map(spliced_df)


def delta_distance_map(
    distances_a: pd.DataFrame, distances_b: pd.DataFrame
) -> alt.Chart:
    df = _merge_pairwise_distances(distances_a, distances_b).assign(
        delta_distance=lambda x: x.distance_a - x.distance_b
    )
    range_max = df.delta_distance.abs().max()
    return (
        alt.Chart(df)
        .mark_rect()
        .encode(
            x=alt.X('atom_id_1', title='atom ID', sort=None),
            y=alt.Y('atom_id_2', title='atom ID', sort=None),
            color=alt.Color(
                'delta_distance',
                scale=alt.Scale(scheme='redblue', domain=[-range_max, range_max]),
            ),
            tooltip=[
                alt.Tooltip('atom_id_1', title='Atom #1'),
                alt.Tooltip('atom_id_2', title='Atom #2'),
                alt.Tooltip('distance_a', title='Apo distance (\u212B)', format='.1f'),
                alt.Tooltip(
                    'distance_b', title='Bound distance (\u212B)', format='.1f'
                ),
                alt.Tooltip(
                    'delta_distance', title='\u0394Distance (\u212B)', format='.1f'
                ),
            ],
        )
        .properties(
            width=600,
            height=600,
        )
    )


def distance_scatter(
    distances_a: pd.DataFrame, distances_b: pd.DataFrame, noe_threshold: float
) -> alt.Chart:
    df = (
        _merge_pairwise_distances(distances_a, distances_b)
        .loc[lambda x: (x.distance_a < noe_threshold) | (x.distance_b < noe_threshold)]
        .assign(delta_distance=lambda x: x.distance_a - x.distance_b)
    )
    range_max = df.delta_distance.abs().max()
    return (
        alt.Chart(df)
        .mark_circle(size=100)
        .encode(
            x=alt.X('distance_a', title='Distance in apo conformation'),
            y=alt.Y('distance_b', title='Distance in bound conformation'),
            color=alt.Color(
                'delta_distance',
                title='\u0394Distance (\u212B)',
                scale=alt.Scale(scheme='redblue', domain=[-range_max, range_max]),
            ),
            opacity=alt.value(0.5),
            tooltip=[
                'atom_id_1',
                'atom_id_2',
                alt.Tooltip('distance_a', format='.1f'),
                alt.Tooltip('distance_b', format='.1f'),
                alt.Tooltip('delta_distance', format='.1f'),
            ],
        )
        .properties(
            width=600,
            height=600,
        )
    )
