import altair as alt
import pandas as pd

from smoltools.rate_my_plate.main import _get_thresholds, _get_time_at_threshold

PLATE_ORDER = [f'{row}{column}' for row in 'ABCDEFGH' for column in range(1, 13)]


def _add_thresholds(
    group: pd.DataFrame, lower_percent: float, upper_percent: float
) -> pd.DataFrame:
    lower_threshold, upper_threshold = _get_thresholds(
        group, lower_percent=lower_percent, upper_percent=upper_percent
    )
    start_time = _get_time_at_threshold(group, lower_threshold)
    end_time = _get_time_at_threshold(group, upper_threshold)
    return group.assign(
        lower_threshold=lower_threshold,
        upper_threshold=upper_threshold,
        start_time=start_time,
        end_time=end_time,
    )


def consumption_curve(
    df: pd.DataFrame, lower_percent: float, upper_percent: float
) -> alt.Chart:
    df_with_thresholds = df.groupby('well', as_index=False).apply(
        _add_thresholds, lower_percent, upper_percent
    )

    scatter = (
        alt.Chart()
        .mark_circle(size=60)
        .encode(
            x=alt.X('time', title='Time (minutes)'),
            y=alt.Y('nadh_consumed', title='NADH consumed'),
            opacity=alt.condition(
                (alt.datum.time >= alt.datum.start_time)
                & (alt.datum.time <= alt.datum.end_time),
                alt.value(0.8),
                alt.value(0.2),
            ),
        )
    ).properties(
        height=100,
        width=150,
    )

    upper_rule = (
        alt.Chart()
        .mark_rule()
        .encode(
            x=alt.X('end_time', title='Time (minutes)'),
            # y=alt.Y('upper_threshold', title='NADH consumed'),
        )
    )
    lower_rule = (
        alt.Chart()
        .mark_rule()
        .encode(
            x=alt.X('start_time', title='Time (minutes)'),
            # y=alt.Y('lower_threshold', title='NADH consumed'),
        )
    )
    return alt.layer(scatter, upper_rule, lower_rule, data=df_with_thresholds).facet(
        facet=alt.Facet('well:N', title='well', sort=PLATE_ORDER),
        columns=6,
    )


def kinetics_curves(df: pd.DataFrame) -> alt.Chart:
    return (
        alt.Chart(df)
        .mark_circle(size=60)
        .encode(
            x=alt.X('column', title='Column'),
            y=alt.Y('rate', title='Rate of NADH consumption / uM protein'),
            tooltip=[
                alt.Tooltip('well', title='Well'),
                alt.Tooltip('rate', title='Rate', format='.2f'),
            ],
            facet=alt.Facet('row', title='Row', columns=4),
        )
        .properties(
            height=200,
            width=250,
        )
    )
