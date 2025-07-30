"""
Functions for systematic imputation / assumption logging.
For now returns df unchanged but logs missing fields.
"""

import pandas as pd, logging

logging.basicConfig(level=logging.INFO, format="IMPUTE: %(message)s")


def impute(df: pd.DataFrame, dataset: str) -> pd.DataFrame:
    # Example: flag missing indication
    missing = df["Indication"].isna().sum()
    if missing:
        logging.info(f"{dataset}: {missing} rows missing 'Indication' – left blank.")
    return df
