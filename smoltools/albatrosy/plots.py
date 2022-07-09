"""Collection of functions for generating plots for TROSY signal."""
import altair as alt
import numpy as np
import pandas as pd


def _distance_map_base(df: pd.DataFrame) -> alt.Chart:
    """Common distance map components."""
    SIZE = 600
    return (
        alt.Chart(df)
        .mark_rect()
        .encode(
            x=alt.X(
                'id_1',
                title='Atom ID [A]',
                sort=None,
                axis=alt.Axis(labels=False, ticks=False),
            ),
            y=alt.Y(
                'id_2',
                title='Atom ID [B]',
                sort=None,
                axis=alt.Axis(labels=False, ticks=False),
            ),
        )
        .properties(
            width=SIZE,
            height=SIZE,
        )
    )


def distance_map(df: pd.DataFrame) -> alt.Chart:
    """Heatmap of pairwise distance between each labelled atom.

    Parameters:
    -----------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.

    Returns:
    --------
    Chart: Altair chart object.
    """
    return _distance_map_base(df).encode(
        color=alt.Color('distance', title='Distance (\u212B)'),
        tooltip=[
            alt.Tooltip('id_1', title='Atom #1'),
            alt.Tooltip('id_2', title='Atom #2'),
            alt.Tooltip('distance', title='Distance (\u212B)', format='.1f'),
        ],
    )


def binned_distance_map(df: pd.DataFrame, bin_size: int) -> alt.Chart:
    """Heatmap of pairwise distance between each labelled atom. Distances are binned
    to reduce visual clutter.

    Parameters:
    -----------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.
    bin_size (int): Bin size for binning distances.

    Returns:
    --------
    Chart: Altair chart object.
    """
    return _distance_map_base(df).encode(
        color=alt.Color(
            'distance', title='Distance (\u212B)', bin=alt.Bin(step=bin_size)
        ),
        tooltip=[
            alt.Tooltip('id_1', title='Atom #1'),
            alt.Tooltip('id_2', title='Atom #2'),
            alt.Tooltip('distance', title='Distance (\u212B)', format='.1f'),
        ],
    )


def _add_noe_bins(df: pd.DataFrame) -> pd.DataFrame:
    """Add column converting distance into relative NOE strength."""
    return df.assign(
        noe_strength=lambda x: pd.cut(
            x.distance,
            bins=[0, 5, 8, 10, np.inf],
            include_lowest=True,
            labels=['strong', 'medium', 'weak', 'none'],
        )
    )


def noe_map(df: pd.DataFrame) -> alt.Chart:
    """Heatmap of expected NOE between each labelled atom.

    Parameters:
    -----------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom pair
        and the distance (in angstroms) between each pair.

    Returns:
    --------
    Chart: Altair chart object.
    """
    return _distance_map_base(df.pipe(_add_noe_bins)).encode(
        color=alt.Color(
            'noe_strength',
            title='NOE',
            scale=alt.Scale(
                domain=['strong', 'medium', 'weak', 'none'],
                scheme='blues',
                reverse=True,
            ),
        ),
        tooltip=[
            alt.Tooltip('id_1', title='Atom #1'),
            alt.Tooltip('id_2', title='Atom #2'),
            alt.Tooltip('distance', title='Distance (\u212B)', format='.1f'),
            alt.Tooltip('noe_strength', title='NOE'),
        ],
    )


def delta_distance_map(df: pd.DataFrame) -> alt.Chart:
    """Heatmap of pairwise distance between each labelled atom.

    Parameters:
    -----------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom
        pair and the distance (in angstroms) between each pair in each of the two
        conformations, as well as the difference in pairwise distance between the
        conformations.

    Returns:
    --------
    Chart: Altair chart object.
    """
    range_max = df.delta_distance.abs().max()

    return _distance_map_base(df).encode(
        color=alt.Color(
            'delta_distance',
            title='\u0394Distance (\u212B)',
            scale=alt.Scale(scheme='redblue', domain=[-range_max, range_max]),
        ),
        tooltip=[
            alt.Tooltip('id_1', title='Atom #1'),
            alt.Tooltip('id_2', title='Atom #2'),
            alt.Tooltip('distance_a', title='Conformation A (\u212B)', format='.1f'),
            alt.Tooltip('distance_b', title='Conformation B (\u212B)', format='.1f'),
            alt.Tooltip(
                'delta_distance', title='\u0394Distance (\u212B)', format='.1f'
            ),
        ],
    )


def distance_scatter(df: pd.DataFrame, noe_threshold: float) -> alt.Chart:
    """Scatter plot of pairwise distance between each labelled atom in one conformation
    versus the other.

    Parameters:
    -----------
    DataFrame: Dataframe with the atom IDs (residue number, carbon ID) of each atom
        pair and the distance (in angstroms) between each pair in each of the two
        conformations, as well as the difference in pairwise distance between the
        conformations.

    Returns:
    --------
    Chart: Altair chart object.
    """
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
                alt.Tooltip('id_1', title='Atom #1'),
                alt.Tooltip('id_2', title='Atom #2'),
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
