import altair as alt
import pandas as pd

from smoltools.rate_my_plate.main import _get_thresholds

PLATE_ORDER = [f'{row}{column}' for row in 'ABCDEFGH' for column in range(1, 13)]


def _add_thresholds(
    df: pd.DataFrame, lower_percent: float, upper_percent: float
) -> pd.DataFrame:
    lower_threshold, upper_threshold = _get_thresholds(
        df, lower_percent=lower_percent, upper_percent=upper_percent
    )
    return df.assign(
        lower_threshold=lower_threshold,
        upper_threshold=upper_threshold,
    )


def consumption_curve(
    df: pd.DataFrame, lower_percent: float, upper_percent: float
) -> alt.Chart:

    lower_threshold, upper_threshold = _get_thresholds(df, lower_percent, upper_percent)

    scatter = (
        alt.Chart()
        .mark_circle(size=60)
        .encode(
            x=alt.X('time', title='Time (minutes)'),
            y=alt.Y('nadh_consumed', title='NADH consumed'),
            opacity=alt.condition(
                (alt.datum.nadh_consumed >= lower_threshold)
                & (alt.datum.nadh_consumed <= upper_threshold),
                alt.value(0.8),
                alt.value(0.2),
            ),
        )
    ).properties(
        height=100,
        width=150,
    )

    upper_rule = (
        alt.Chart(pd.DataFrame({'upper_threshold': [upper_threshold]}))
        .mark_rule()
        .encode(y=alt.Y('upper_threshold', title='NADH consumed'))
    )
    lower_rule = (
        alt.Chart(pd.DataFrame({'lower_threshold': [lower_threshold]}))
        .mark_rule()
        .encode(y=alt.Y('lower_threshold', title='NADH consumed'))
    )
    return alt.layer(scatter, upper_rule, lower_rule, data=df).facet(
        facet=alt.Facet('well:N', title='well', sort=PLATE_ORDER),
        columns=6,
    )


def kinetics_curves(df: pd.DataFrame) -> alt.Chart:
    return (
        alt.Chart(df)
        .mark_circle(size=60)
        .encode(
            x=alt.X('column', title='Plate column'),
            y=alt.Y('rate', title='Rate of NADH consumption'),
            facet=alt.Facet('row', columns=4),
        )
        .properties(
            height=150,
            width=200,
        )
    )
