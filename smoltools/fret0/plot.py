"""Collection of functions for generating plots for finding candidate residues for
FRET."""

import altair as alt
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from fret0.utils import convert_to_path
from fret0.efficiency import _generate_r0_curve

LIGHT_GREY = '#969aa8'
RED = '#d62728'
BLUE = '#1f77b4'


def _pivot_df(df: pd.DataFrame, pivot_value: str) -> pd.DataFrame:
    """Pivots long form DataFrame to wide form for heatmap plotting."""
    return df.pivot(
        columns='residue_number_1', index='residue_number_2', values=pivot_value
    )


def _delta_heatmap(df: pd.DataFrame, title: str, ax: plt.Axes = None) -> None:
    """Generate a heatmap encoding the difference in a value between conformations"""
    # TODO: replace w/ altair
    with sns.axes_style('ticks'):
        if ax is None:
            plt.figure(figsize=(20, 16))
            sns.heatmap(df)
            plt.xlabel('Residue Number', fontsize=12)
            plt.ylabel('Residue Number', fontsize=12)
            plt.title(title, fontsize=16)
        else:
            sns.heatmap(df, ax=ax)
            ax.set_xlabel('Residue Number', fontsize=9)
            ax.set_ylabel('Residue Number', fontsize=9)
            ax.set_title(title, fontsize=11)


def distance_map(df: pd.DataFrame, ax: plt.Axes = None) -> None:
    """
    Generates a heatmap encoding the difference in pairwise residue distances between
    two conformations. Optionally exports the result to a .png if a path is provided.

    Args:
        df (DataFrame): dataframe with difference in pairwise distances between
            conformations (long form).
        export_path (Path): Optional export path for .png of the heatmap.

    Returns:
        None
    """
    df = _pivot_df(df, 'delta_distance')
    title = 'Pairwise Δdistance (Å) between conformations'
    _delta_heatmap(df, title=title, ax=ax)


def e_fret_map(df: pd.DataFrame, ax: plt.Axes = None) -> None:
    """
    Generates a heatmap encoding the difference in residue pair FRET efficiencies
    between two conformations. Optionally exports the result to a .png if a path is
    provided.

    Args:
        df (DataFrame): dataframe with difference in FRET effciencies between residue
            pairs in two conformations (long form).
        export_path (Path): Optional export path for .png of the heatmap.

    Returns:
        None
    """
    r0 = df.r0[0]
    df = _pivot_df(df, 'delta_E_fret')
    title = f'Pairwise ΔE_fret (@r0={r0}) between conformations'
    _delta_heatmap(df, title=title, ax=ax)


def r0_curves(distance_a: float, distance_b: float) -> alt.Chart:
    """
    Generates an interactive plot to visualize the FRET efficiencies for two residue
    pair differences as a function of R0.

    Args:
        distance_a (float): distance between FRET donor and acceptor in conformation_a
        distance_b (float): distance between FRET donor and acceptor in conformation_b

    Returns:
        Chart: interactive Altair chart.
    """
    e_fret_by_distance, e_fret_delta = _generate_r0_curve(distance_a, distance_b)
    nearest = alt.selection(
        type='single', nearest=True, on='mouseover', fields=['r0'], empty='none'
    )

    line = (
        alt.Chart(
            e_fret_by_distance,
        )
        .mark_line(interpolate='basis')
        .encode(
            x=alt.X('r0', title='R0 (\u212B)'),
            y=alt.Y('e_fret', title='E_fret'),
            color=alt.Color(
                'distance',
                scale=alt.Scale(domain=['A', 'B'], range=[RED, BLUE]),
                legend=alt.Legend(orient='bottom-right'),
            ),
        )
    )

    delta = (
        alt.Chart(e_fret_delta)
        .mark_area(interpolate='basis')
        .encode(
            x='r0',
            y='delta',
            opacity=alt.value(0.3),
            color=alt.value(LIGHT_GREY),
        )
    )

    selectors = (
        alt.Chart(e_fret_by_distance)
        .mark_point()
        .encode(
            x='r0',
            opacity=alt.value(0),
        )
        .add_selection(nearest)
    )

    points = line.mark_circle(size=50).encode(
        opacity=alt.condition(nearest, alt.value(1), alt.value(0))
    )

    text = (
        line.mark_text(align='left', dx=10, dy=10, fontWeight='bold')
        .encode(text=alt.condition(nearest, 'label:O', alt.value('')))
        .transform_calculate(label='format(datum.e_fret,".1%")')
    )

    rules = (
        alt.Chart(e_fret_by_distance)
        .mark_rule(color='gray')
        .encode(
            x='r0',
        )
        .transform_filter(nearest)
    )

    return alt.layer(delta, line, selectors, points, rules, text)


def _export_figure(export_path) -> None:
    export_path = convert_to_path(export_path)
    plt.savefig(export_path.with_suffix('.png'), dpi=300)
