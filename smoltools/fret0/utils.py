import pandas as pd


def splice_e_fret_tables(df: pd.DataFrame) -> pd.DataFrame:
    df.loc[lambda x: x.atom_id_1 <= x.atom_id_2, 'E_fret'] = df['E_fret_a']
    df.loc[lambda x: x.atom_id_1 > x.atom_id_2, 'E_fret'] = df['E_fret_b']
    return df[['atom_id_1', 'atom_id_2', 'E_fret', 'delta_E_fret']]
