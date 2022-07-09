import pandas as pd


def extract_residue_number(s: pd.Series) -> pd.Series:
    return s.str.partition('-')[0].str[3:].astype(int)


def splice_conformation_tables(df_a: pd.DataFrame, df_b: pd.DataFrame) -> pd.DataFrame:
    """Splice distance tables for two conformations together.

    Parameters:
    -----------
    df_a (DataFrame): Dataframe with the atom IDs (residue number, carbon ID) of each
        atom pair and the distance (in angstroms) between each pair.
    df_b (DataFrame): Dataframe with the atom IDs (residue number, carbon ID) of each
        atom pair and the distance (in angstroms) between each pair.

    Returns:
    --------
    DataFrame: DataFrame with the lower triangle of the DataFrame containing values
        from the first conformation and the upper triangle of the DataFrame containing
        values from the second conformation.
    """
    return pd.concat(
        [
            df_a.loc[
                lambda x: extract_residue_number(x.id_1)
                <= extract_residue_number(x.id_2)
            ],
            df_b.loc[
                lambda x: extract_residue_number(x.id_1)
                > extract_residue_number(x.id_2)
            ],
        ]
    ).sort_values(['id_1', 'id_2'], key=extract_residue_number)
