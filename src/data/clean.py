"""
Full cleaning pipeline:
* standardise dates
* normalise drug names / indications
* enforce column schema
* split combination drugs into individual APIs
"""

from __future__ import annotations
import pandas as pd, re, unicodedata
import argparse
from .load import load_raw
from .impute import impute
from ..config import CDSCO_CLEAN, FDA_CLEAN, CDSCO_EXPLODED, FDA_VOCAB_PATH, PROC
from ..utils.api_vocab import build_and_save_vocabulary, load_vocabulary
from ..utils.api_splitter import split_apis

PROC.mkdir(exist_ok=True, parents=True)

DATE_PAT = re.compile(r"(\d{1,2})[/-](\d{1,2})[/-](\d{2,4})")


def _norm_date(s: str) -> str | pd.NaT:
    """Normalize date strings to a consistent format (MM/DD/YYYY)."""
    if pd.isna(s) or not str(s).strip():
        return pd.NaT
    s = str(s).strip()
    if re.match(r"^\d{4}$", s):  # year only
        return f"01/01/{s}"
    if re.match(r"^\d{2}/\d{4}$", s):  # mm/yyyy
        m, y = s.split("/")
        return f"{m.zfill(2)}/01/{y}"
    m = DATE_PAT.search(s.replace(".", "/").replace("-", "/"))
    if m:
        d, mth, yr = m.groups()
        yr = yr.zfill(4) if len(yr) == 4 else f"19{yr}" if int(yr) > 30 else f"20{yr}"
        return f"{mth.zfill(2)}/{d.zfill(2)}/{yr}"
    return pd.NaT


def _norm_text(x: str) -> str:
    """Normalize text by removing extra spaces and accents."""
    if pd.isna(x):
        return ""
    x = unicodedata.normalize("NFKD", x)
    x = re.sub(r"\s+", " ", x).strip()
    return x


def _strip_forms(drug: str) -> str:
    """Remove dosage forms and strengths from drug names."""
    drug = _norm_text(drug)
    drug = re.sub(
        r"\b(tablet|capsule|injection|cream|ointment|spray|solution|gel|drops?|suspension|eye|ear|nasal|intranasal|oral|iv|im|vial|ampoule|sachet|mg|mcg|g|%|w\/v|w\/w|v\/v)\b.*",
        "",
        drug,
        flags=re.I,
    )
    drug = re.sub(r"[\d.,]+(mg|mcg|g|%)+.*", "", drug, flags=re.I)
    drug = drug.strip(",;- ").title()
    return drug


def _explode_combination_drugs(df: pd.DataFrame) -> pd.DataFrame:
    """
    Explode combination drugs into individual API components.
    
    Args:
        df: DataFrame with cleaned drug data
        
    Returns:
        DataFrame with one row per API component, preserving metadata
    """
    # Load or build FDA vocabulary
    vocab = load_vocabulary(FDA_VOCAB_PATH)
    if not vocab:
        print("Building FDA API vocabulary...")
        vocab = build_and_save_vocabulary()
    
    if not vocab:
        print("⚠️ Warning: Empty FDA vocabulary, can't split combination drugs")
        return df
    
    # Apply API splitting to each drug name
    print("Splitting combination drugs into individual APIs...")
    exploded_rows = []
    
    for _, row in df.iterrows():
        drug_name = row["Drug Name"]
        apis = split_apis(drug_name, vocab)
        
        if not apis:
            # If no APIs found, keep the original drug name as a single API
            exploded_rows.append(row)
        else:
            # Create a new row for each API
            for api in apis:
                new_row = row.copy()
                new_row["Drug Name"] = api
                new_row["Original Drug Name"] = drug_name  # Store the original combination
                new_row["Is Combination"] = len(apis) > 1  # Flag if this was part of a combination
                exploded_rows.append(new_row)
    
    exploded_df = pd.DataFrame(exploded_rows)
    
    # Add normalized drug names again after splitting
    if "drug_norm" in exploded_df.columns:
        from ..utils.text import normalize
        exploded_df["drug_norm"] = exploded_df["Drug Name"].apply(normalize)
    
    return exploded_df


def clean(explode_combinations: bool = False):
    """
    Load, clean, and standardize CDSCO and FDA datasets, then save cleaned data.
    
    Args:
        explode_combinations: If True, explode combination drugs into individual APIs
    """
    try:
        cdsco, fda = load_raw()
    except FileNotFoundError as e:
        print(f"Error loading raw data: {e}")
        return

    # Standardize CDSCO columns
    cdsco = cdsco.rename(
        columns={
            "Drug Name": "Drug Name",
            "Indication": "Indication",
            "Date of Approval": "Date of Approval",
            "Strength": "Strength",
        }
    )
    cdsco["Date of Approval"] = cdsco["Date of Approval"].apply(_norm_date)
    cdsco["Drug Name"] = cdsco["Drug Name"].apply(_strip_forms)
    cdsco["Indication"] = cdsco["Indication"].apply(_norm_text)

    # Log missing fields and impute if necessary for CDSCO
    cdsco = impute(cdsco, "CDSCO")

    # Standardize FDA columns to match CDSCO schema
    fda = fda.rename(
        columns={
            "Generic Name": "Drug Name",
            "Approved Labeled Indication": "Indication",
            "Marketing Approval Date": "Date of Approval",
            "Orphan Designation": "Orphan",
            "Sponsor Company": "Sponsor",
            "Sponsor Country": "Country",
        }
    )
    fda["Date of Approval"] = fda["Date of Approval"].apply(_norm_date)
    fda["Drug Name"] = fda["Drug Name"].apply(_strip_forms)
    fda["Indication"] = fda["Indication"].apply(_norm_text)
    if "Orphan" in fda:
        fda["Orphan"] = (
            fda["Orphan"]
            .astype(str)
            .str.lower()
            .str.contains("yes|designat|approved")
            .astype(int)
        )

    # Log missing fields and impute if necessary for FDA
    fda = impute(fda, "FDA")

    # Drop duplicate drug entries to minimize duplicates after cleaning
    cdsco = cdsco.drop_duplicates(subset=["Drug Name"])
    fda = fda.drop_duplicates(subset=["Drug Name"])

    # Enrich with RxNorm IDs where possible
    try:
        from ..utils.rxnorm import name_to_rxcui
        from ..utils.text import normalize as txt_normalize
        print("🔍 Enriching with RxNorm IDs...")
        cdsco["RxCUI"] = cdsco["Drug Name"].apply(lambda x: name_to_rxcui(txt_normalize(x)))
        fda["RxCUI"] = fda["Drug Name"].apply(lambda x: name_to_rxcui(txt_normalize(x)))
    except ImportError:
        print("⚠️ RxNorm enrichment not available; continuing without RxCUI.")

    # Save standard cleaned data (with RxCUI)
    cdsco.to_parquet(CDSCO_CLEAN, index=False)
    fda.to_parquet(FDA_CLEAN, index=False)
    print("✅ cleaned data written to", PROC)
    
    # Process combination drugs if requested
    if explode_combinations:
        # First save FDA vocabulary for API matching
        from ..utils.api_vocab import build_and_save_vocabulary
        build_and_save_vocabulary()
        
        # Now explode CDSCO combination drugs
        cdsco_exploded = _explode_combination_drugs(cdsco)
        cdsco_exploded.to_parquet(CDSCO_EXPLODED, index=False)
        print(f"✅ Exploded {len(cdsco)} CDSCO entries into {len(cdsco_exploded)} individual APIs")
        print(f"  Saved to {CDSCO_EXPLODED}")


def main():
    """Parse command-line arguments and run the data cleaning process."""
    parser = argparse.ArgumentParser(
        description="Clean and standardize drug data from CDSCO and FDA"
    )
    parser.add_argument(
        "--explode-combinations",
        action="store_true",
        help="Explode combination drugs into individual APIs",
    )
    args = parser.parse_args()
    clean(explode_combinations=args.explode_combinations)


if __name__ == "__main__":
    main()
