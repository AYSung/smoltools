import altair as alt
import numpy as np
import pandas as pd


def _distance_map_base(df: pd.DataFrame) -> alt.Chart:
    SIZE = 600
    return (
        alt.Chart(df)
        .mark_rect()
        .encode(
            x=alt.X('atom_id_1', title='Atom ID', sort=None),
            y=alt.Y('atom_id_2', title='Atom ID', sort=None),
        )
        .properties(
            width=SIZE,
            height=SIZE,
        )
    )


def distance_map(df: pd.DataFrame) -> alt.Chart:
    return _distance_map_base(df).encode(
        color=alt.Color('distance', title='Distance (\u212B)'),
        tooltip=[
            alt.Tooltip('atom_id_1', title='Atom #1'),
            alt.Tooltip('atom_id_2', title='Atom #2'),
            alt.Tooltip('distance', title='Distance (\u212B)', format='.1f'),
        ],
    )


def binned_distance_map(df: pd.DataFrame, bin_size: int) -> alt.Chart:
    return _distance_map_base(df).encode(
        color=alt.Color(
            'distance', title='Distance (\u212B)', bin=alt.Bin(step=bin_size)
        ),
        tooltip=[
            alt.Tooltip('atom_id_1', title='Atom #1'),
            alt.Tooltip('atom_id_2', title='Atom #2'),
            alt.Tooltip('distance', title='Distance (\u212B)', format='.1f'),
        ],
    )


def _add_noe_bins(df: pd.DataFrame) -> pd.DataFrame:
    return df.assign(
        noe_strength=lambda x: pd.cut(
            x.distance,
            bins=[0, 2.5, 3.5, 5, 15, np.inf],
            include_lowest=True,
            labels=['strong', 'medium', 'weak', 'very weak', 'none'],
        )
    )


def noe_map(df: pd.DataFrame) -> alt.Chart:
    return _distance_map_base(df.pipe(_add_noe_bins)).encode(
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


def delta_distance_map(df: pd.DataFrame) -> alt.Chart:
    range_max = df.delta_distance.abs().max()

    return _distance_map_base(df).encode(
        color=alt.Color(
            'delta_distance',
            title='\u0394Distance (\u212B)',
            scale=alt.Scale(scheme='redblue', domain=[-range_max, range_max]),
        ),
        tooltip=[
            alt.Tooltip('atom_id_1', title='Atom #1'),
            alt.Tooltip('atom_id_2', title='Atom #2'),
            alt.Tooltip('distance_a', title='Conformation A (\u212B)', format='.1f'),
            alt.Tooltip('distance_b', title='Conformation B (\u212B)', format='.1f'),
            alt.Tooltip(
                'delta_distance', title='\u0394Distance (\u212B)', format='.1f'
            ),
        ],
    )


def distance_scatter(df: pd.DataFrame, noe_threshold: float) -> alt.Chart:
    range_max = df.delta_distance.abs().max()

    return (
        alt.Chart(
            df.loc[
                lambda x: (x.distance_a < noe_threshold)
                | (x.distance_b < noe_threshold)
            ]
        )
        .mark_circle(size=100)
        .encode(
            x=alt.X('distance_a', title='Distance in Conformation A'),
            y=alt.Y('distance_b', title='Distance in Conformation B'),
            color=alt.Color(
                'delta_distance',
                title='\u0394Distance (\u212B)',
                scale=alt.Scale(scheme='redblue', domain=[-range_max, range_max]),
            ),
            opacity=alt.value(0.5),
            tooltip=[
                alt.Tooltip('atom_id_1', title='Atom #1'),
                alt.Tooltip('atom_id_2', title='Atom #2'),
                alt.Tooltip(
                    'distance_a', title='Conformation A (\u212B)', format='.1f'
                ),
                alt.Tooltip(
                    'distance_b', title='Conformation B (\u212B)', format='.1f'
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
