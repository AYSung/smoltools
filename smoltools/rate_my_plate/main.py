from io import BytesIO
from pathlib import Path

import pandas as pd
from scipy.stats import linregress


def _time_to_minutes(time: str) -> float:
    hours, minutes, seconds = time.split(":")
    return (60 * int(hours)) + int(minutes) + (int(seconds) / 60)


def read_data(path: str) -> pd.DataFrame:
    if not isinstance(path, Path):
        path = Path(path)

    return pd.read_excel(path, skiprows=2).pipe(clean_import)


def read_data_from_bytes(bytes_data: bytes) -> pd.DataFrame:
    return pd.read_excel(BytesIO(bytes_data), skiprows=1).pipe(clean_import)


def clean_import(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.rename(columns={"Kinetic read": "time"})
        .astype({"time": str})
        .pipe(convert_time)
        .pipe(absorbance_to_consumption)
        .pipe(tidy_data)
    )


def convert_time(df: pd.DataFrame) -> pd.DataFrame:
    """convert time columns to fractions of a minute and set as index."""
    return df.assign(time=lambda x: x.time.map(_time_to_minutes)).set_index("time")


def absorbance_to_consumption(df: pd.DataFrame) -> pd.DataFrame:
    """
    convert absorbance values into NADH consumption by normalizing to time zero and
    dividing by extinction coefficient of NADH (0.00622).
    """
    NADH_EXTINCTION_COEFFICIENT = 0.00622

    return (df.iloc[0] - df) / NADH_EXTINCTION_COEFFICIENT


def tidy_data(df: pd.DataFrame) -> pd.DataFrame:
    return df.melt(
        var_name="well", value_name="nadh_consumed", ignore_index=False
    ).reset_index()


def _get_thresholds(
    df: pd.DataFrame, lower_percent, upper_percent
) -> tuple[float, float]:
    max_value = df.nadh_consumed.max()
    lower_threshold = lower_percent * max_value
    upper_threshold = upper_percent * max_value
    return lower_threshold, upper_threshold


def _get_time_at_threshold(df: pd.DataFrame, threshold: float) -> float:
    if threshold > df.nadh_consumed.max():
        return df.time.max()
    elif threshold < df.nadh_consumed.min():
        return 0
    else:
        return df.loc[lambda x: x.nadh_consumed >= threshold].time.min()


def filter_data(
    df: pd.DataFrame, lower_percent: float, upper_percent: float
) -> pd.DataFrame:
    def _filter_well(group: pd.DataFrame) -> pd.DataFrame:
        lower_threshold, upper_threshold = _get_thresholds(
            group, lower_percent, upper_percent
        )
        start_time = _get_time_at_threshold(group, lower_threshold)
        end_time = _get_time_at_threshold(group, upper_threshold)
        return group.loc[lambda x: (x.time >= start_time) & (x.time <= end_time)]

    return df.groupby("well", as_index=False).apply(_filter_well).reset_index(drop=True)


def _estimate_slope(df: pd.DataFrame) -> float:
    slope = linregress(df.time, df.nadh_consumed)[0]
    return slope


def calculate_slopes(df: pd.DataFrame) -> pd.DataFrame:
    """Calculates rate of NADH consumption / ATP production through linear regression."""
    return (
        df.groupby("well", as_index=False)
        .apply(_estimate_slope)
        .rename(columns={None: "rate"})
        .assign(
            row=lambda x: x.well.str[:1], column=lambda x: x.well.str[1:].astype(int)
        )
        .sort_values(["row", "column"])
    )


def normalize_to_protein_concentration(
    df: pd.DataFrame, concentration: float
) -> pd.DataFrame:
    """Normalizes NADH consumption / ATP production rate to provided protein concentration (in uM)."""
    return df.assign(rate=lambda x: x.rate / concentration)


def convert_to_wide(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.pivot(index="column", columns="row", values="rate")
        .reset_index()
        .rename_axis(columns=None)
    )


def rate_plate(
    df: pd.DataFrame,
    lower_percent: float,
    upper_percent: float,
    concentration: float = 1,
) -> pd.DataFrame:
    return (
        df.pipe(filter_data, lower_percent=lower_percent, upper_percent=upper_percent)
        .pipe(calculate_slopes)
        .pipe(normalize_to_protein_concentration, concentration=concentration)
    )
