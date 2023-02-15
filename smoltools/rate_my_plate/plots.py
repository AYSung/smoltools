import altair as alt
import pandas as pd

from smoltools.rate_my_plate.main import _get_thresholds


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

    df_with_thresholds = _add_thresholds(df, lower_percent, upper_percent)
    scatter = (
        alt.Chart()
        .mark_circle(size=60)
        .encode(
            x=alt.X('time', title='Time (minutes)'),
            y=alt.Y('nadh_consumed', title='NADH consumed'),
            facet=alt.Facet('well', columns=12),
        )
    )
    upper_rule = (
        alt.Chart()
        .mark_rule()
        .encode(y=alt.Y('upper_threshold', title='NADH consumed'))
    )
    lower_rule = (
        alt.Chart()
        .mark_rule()
        .encode(y=alt.Y('lower_threshold', title='NADH consumed'))
    )
    return alt.layer(
        scatter, upper_rule, lower_rule, data=df_with_thresholds
    ).properties(
        height=100,
        width=150,
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
