import pandas as pd


def extract_residue_number(s: pd.Series) -> pd.Series:
    return s.str[3:].astype(int)


def lower_triangle(df: pd.DataFrame) -> pd.Series:
    return extract_residue_number(df.id_1) < extract_residue_number(df.id_2)


def sort_table(df: pd.DataFrame) -> pd.DataFrame:
    return df.sort_values(['id_1', 'id_2'], key=extract_residue_number)
