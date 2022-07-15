"""Collection of functions for generating plots for finding candidate residues for
FRET."""

import altair as alt
import pandas as pd

from smoltools.fret0.efficiency import generate_r0_curve
import smoltools.resources.colors as colors
from smoltools.fret0.utils import lower_triangle, sort_table


def _get_size(n_residues: int) -> int:
    MAX_SIZE = 600
    return min(MAX_SIZE, n_residues * 10)


def _distance_map_base(df: pd.DataFrame) -> alt.Chart:
    """Common distance map components."""
    df = sort_table(df)

    n_residues = df.id_1.nunique()
    size = _get_size(n_residues)

    axis_config = {'sort': None, 'axis': alt.Axis(labels=False, ticks=False)}

    return (
        alt.Chart(df)
        .mark_rect()
        .encode(
            x=alt.X('id_1:N', title='Residue #', **axis_config),
            y=alt.Y('id_2:N', title='Residue #', **axis_config),
        )
        .properties(width=size, height=size)
    )


def delta_distance_map(df: pd.DataFrame, cutoff: float = 5) -> alt.Chart:
    """Heatmap of pairwise distance between each alpha carbon between two conformations.

    Parameters:
    -----------
    DataFrame: Dataframe with the atom IDs (residue number) of each alpha carbon pair
        and the distance (in angstroms) between each pair.

    Returns:
    --------
    Chart: Altair chart object.
    """
    df = df.loc[lambda x: lower_triangle(x) & (x.delta_distance.abs() > cutoff)]

    range_max = df.delta_distance.abs().max()

    return _distance_map_base(df).encode(
        color=alt.Color(
            'delta_distance',
            title='\u0394Distance (\u212B)',
            scale=alt.Scale(domain=[-range_max, range_max], scheme='redblue'),
        ),
        tooltip=[
            alt.Tooltip('id_1', title='Residue #1'),
            alt.Tooltip('id_2', title='Residue #2'),
            alt.Tooltip('distance_a', title='Conformation A (\u212B)', format='.1f'),
            alt.Tooltip('distance_b', title='Conformation B (\u212B)', format='.1f'),
            alt.Tooltip(
                'delta_distance', title='\u0394Distance (\u212B)', format='.1f'
            ),
        ],
    )


def delta_e_fret_map(df: pd.DataFrame, cutoff: float = 0.1) -> alt.Chart:
    """Heatmap of the difference in E_fret between each alpha carbon between two
    conformations.

    Parameters:
    -----------
    DataFrame: Dataframe with the atom IDs (residue number) of each atom pair and the
        E_fret between each pair in each of the two conformations, as well as the
        difference in the E_fret of each pair between the conformations.

    Returns:
    --------
    Chart: Altair chart object.
    """
    range_max = df.delta_E_fret.abs().max()

    return _distance_map_base(
        df.loc[lambda x: lower_triangle(x) & (x.delta_E_fret.abs() > cutoff)]
    ).encode(
        color=alt.Color(
            'delta_E_fret',
            title='\u0394E_fret',
            scale=alt.Scale(domain=[-range_max, range_max], scheme='redblue'),
        ),
        tooltip=[
            alt.Tooltip('id_1', title='Residue #1'),
            alt.Tooltip('id_2', title='Residue #2'),
            alt.Tooltip('E_fret_a', title='Conformation A', format='.2f'),
            alt.Tooltip('E_fret_b', title='Conformation B', format='.2f'),
            alt.Tooltip('delta_E_fret', title='\u0394E_fret', format='.2f'),
        ],
    )


def e_fret_scatter(df: pd.DataFrame, cutoff: float = 0.2) -> alt.Chart:
    """Scatter plot of pairwise E_fret between each alpha carbon in one conformation
    versus the other.

    Parameters:
    -----------
    DataFrame: Dataframe with the atom IDs (residue number) of each atom pair and the
        E_fret between each pair in each of the two conformations, as well as the
        difference in the E_fret of each pair between the conformations.

    Returns:
    --------
    Chart: Altair chart object.
    """
    range_max = df.delta_E_fret.abs().max()

    return (
        alt.Chart(df.loc[lambda x: x.delta_E_fret.abs() > cutoff])
        .mark_circle(size=100)
        .encode(
            x=alt.X('E_fret_a', title='E_fret in A'),
            y=alt.Y('E_fret_b', title='E_fret in B'),
            color=alt.Color(
                'delta_E_fret',
                title='\u0394E_fret',
                scale=alt.Scale(domain=[-range_max, range_max], scheme='redblue'),
            ),
            opacity=alt.value(0.4),
            tooltip=[
                alt.Tooltip('id_1', title='Residue #1'),
                alt.Tooltip('id_2', title='Residue #2'),
            ],
        )
        .properties(width=600, height=600)
    )


def r0_curves(distance_a: float, distance_b: float) -> alt.Chart:
    """
    Generates an interactive plot to visualize the FRET efficiencies for two residue
    pair differences as a function of R0.

    Parameters:
    -----------
    distance_a (float): distance between FRET donor and acceptor in conformation_a
    distance_b (float): distance between FRET donor and acceptor in conformation_b

    Returns:
    --------
    Chart: Altair chart object.
    """
    e_fret_by_distance, e_fret_delta = generate_r0_curve(distance_a, distance_b)
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
                scale=alt.Scale(domain=['A', 'B'], range=[colors.RED, colors.BLUE]),
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
            color=alt.value(colors.LIGHT_GREY),
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
