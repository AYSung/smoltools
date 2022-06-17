import pandas as pd


def extract_residue_number(s: pd.Series) -> pd.Series:
    return s.str.partition('-')[0].astype(int)


def splice_conformation_tables(df_a: pd.DataFrame, df_b: pd.DataFrame) -> pd.DataFrame:
    return pd.concat(
        [
            df_a.loc[
                lambda x: extract_residue_number(x.atom_id_1)
                <= extract_residue_number(x.atom_id_2)
            ],
            df_b.loc[
                lambda x: extract_residue_number(x.atom_id_1)
                > extract_residue_number(x.atom_id_2)
            ],
        ]
    ).sort_values(['atom_id_1', 'atom_id_2'], key=extract_residue_number)
