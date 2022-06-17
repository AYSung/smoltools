"""Collection of functions for generating plots for finding candidate residues for
FRET."""

import altair as alt
import pandas as pd

from smoltools.fret0.efficiency import generate_r0_curve
import smoltools.calculate.distance as distance
import smoltools.resources.colors as colors


def _distance_map_base(df: pd.DataFrame) -> alt.Chart:
    SIZE = 600
    return (
        alt.Chart(df)
        .mark_rect()
        .encode(
            x=alt.X(
                'atom_id_1:O',
                title='Residue #',
                axis=alt.Axis(labels=False, ticks=False),
            ),
            y=alt.Y(
                'atom_id_2:O',
                title='Residue #',
                axis=alt.Axis(labels=False, ticks=False),
            ),
        )
        .properties(
            width=SIZE,
            height=SIZE,
        )
    )


def delta_distance_map(
    distances_a: pd.DataFrame, distances_b: pd.DataFrame, cutoff: int = 5
) -> alt.Chart:
    df = (
        distance._merge_pairwise_distances(distances_a, distances_b)
        .assign(delta_distance=lambda x: (x.distance_a - x.distance_b))
        .loc[lambda x: (x.atom_id_1 < x.atom_id_2) & (x.delta_distance.abs() > cutoff)]
    )

    range_max = df.delta_distance.abs().max()

    return _distance_map_base(df).encode(
        color=alt.Color(
            'delta_distance',
            title='\u0394Distance (\u212B)',
            scale=alt.Scale(domain=[-range_max, range_max], scheme='redblue'),
        ),
        tooltip=[
            alt.Tooltip('atom_id_1', title='Residue #1'),
            alt.Tooltip('atom_id_2', title='Residue #2'),
            alt.Tooltip('distance_a', title='Conformation A (\u212B)', format='.1f'),
            alt.Tooltip('distance_b', title='Conformation B (\u212B)', format='.1f'),
            alt.Tooltip(
                'delta_distance', title='\u0394Distance (\u212B)', format='.1f'
            ),
        ],
    )


def delta_e_fret_map(df: pd.DataFrame, cutoff: float = 0.1) -> alt.Chart:
    range_max = df.delta_E_fret.abs().max()

    return _distance_map_base(
        df.loc[lambda x: (x.atom_id_1 < x.atom_id_2) & (x.delta_E_fret.abs() > cutoff)]
    ).encode(
        color=alt.Color(
            'delta_E_fret',
            title='\u0394E_fret',
            scale=alt.Scale(domain=[-range_max, range_max], scheme='redblue'),
        ),
        tooltip=[
            alt.Tooltip('atom_id_1', title='Residue #1'),
            alt.Tooltip('atom_id_2', title='Residue #2'),
            alt.Tooltip('E_fret_a', title='Conformation A', format='.2f'),
            alt.Tooltip('E_fret_b', title='Conformation B', format='.2f'),
            alt.Tooltip('delta_E_fret', title='\u0394E_fret', format='.2f'),
        ],
    )


def e_fret_scatter(df: pd.DataFrame, cutoff: float = 0.2) -> alt.Chart:
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
                alt.Tooltip('atom_id_1', title='Residue #1'),
                alt.Tooltip('atom_id_2', title='Residue #2'),
            ],
        )
        .properties(width=600, height=600)
    )


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
