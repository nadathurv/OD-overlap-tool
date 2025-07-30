import argparse
import logging
from pathlib import Path

import pandas as pd
import textdistance
from ..config import CDSCO_CLEAN, FDA_CLEAN, CDSCO_EXPLODED, PROC
from ..utils.text import normalize, jaccard
from rapidfuzz import fuzz
from ..utils.synonyms import load_synonyms

# ----------------------------------------------------------------------------
# Similarity thresholds (adjust for precision/recall tradeoff)
JARO_THRESHOLD = 0.70       # Optimized: was 0.85, now +4.9% recall
JACCARD_THRESHOLD = 0.1     # Minimum Jaccard similarity for candidate filtering
TOKEN_THRESHOLD = 80        # Minimum RapidFuzz token-set ratio for fuzzy match  
RATIO_THRESHOLD = 80        # Minimum Levenshtein ratio for fuzzy match

def is_high_confidence_match(jw, token, ratio, jw_thresh, token_thresh, ratio_thresh):
    """
    Determine if a candidate pair is a strong match by requiring at least two of three
    similarity metrics to meet or exceed their respective thresholds.
    """
    criteria = [
        jw >= jw_thresh,
        token >= token_thresh,
        ratio >= ratio_thresh,
    ]
    return sum(criteria) >= 2
# ----------------------------------------------------------------------------

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

# Ensure output directory exists
PROC.mkdir(exist_ok=True, parents=True)

# Helper for robust fuzzy matching


def jaro(s1, s2):
    """Calculate Jaro-Winkler similarity between two strings."""
    if not isinstance(s1, str) or not isinstance(s2, str):
        return 0.0
    return textdistance.jaro_winkler(s1, s2)


def run(
    threshold: float = JARO_THRESHOLD,
    jaccard_threshold: float = JACCARD_THRESHOLD,
    token_threshold: int = TOKEN_THRESHOLD,
    ratio_threshold: int = RATIO_THRESHOLD,
    out_file: Path | str | None = None,
    use_exploded: bool = False,
):
    """Identify overlapping drugs between CDSCO and FDA datasets using fuzzy matching.
    threshold: Jaro-Winkler similarity threshold.
    jaccard_threshold: Jaccard similarity threshold for candidate filtering.
    token_threshold: Token-set ratio threshold for fuzzy matching.
    ratio_threshold: Levenshtein ratio threshold for fuzzy matching.
    out_file: output CSV file path.
    use_exploded: If True, use the exploded CDSCO dataset with individual APIs."""
    # Determine output path
    if out_file is None:
        out_file = PROC / "overlap.csv"
    out_file = Path(out_file)
    out_file.parent.mkdir(parents=True, exist_ok=True)

    try:
        logging.info("Loading cleaned datasets")
        # Choose CDSCO dataset based on use_exploded flag
        cdsco_path = CDSCO_EXPLODED if use_exploded else CDSCO_CLEAN
        logging.info(f"Using CDSCO dataset: {cdsco_path}")
        
        if not cdsco_path.exists():
            if use_exploded:
                logging.error(f"Exploded CDSCO file not found. Run clean.py with --explode-combinations first.")
                return
            else:
                logging.error(f"CDSCO clean file not found. Run clean.py first.")
                return
                
        c = pd.read_parquet(cdsco_path)
        f = pd.read_parquet(FDA_CLEAN)
        # Initialize match container
        matches = []
        
        # Check for RxCUI-based authoritative matches
        if 'RxCUI' in c.columns and 'RxCUI' in f.columns:
            logging.info("Performing RxNorm ID-based matching")
            id_c = c[c['RxCUI'].notna()]
            id_f = f[f['RxCUI'].notna()]
            id_pairs = id_c.merge(id_f, on='RxCUI', suffixes=('_cdsco','_fda'))
            for _, ip in id_pairs.iterrows():
                matches.append({
                    'CDSCO Drug Name': ip['Drug Name_cdsco'],
                    'FDA Drug Name': ip['Drug Name_fda'],
                    'Similarity Score': 1.0,
                    'Match Type': 'RxNorm',
                    'CDSCO Approval Date': ip.get('Date of Approval_cdsco',''),
                    'FDA Approval Date': ip.get('Date of Approval_fda',''),
                    'CDSCO Indication': ip.get('Indication_cdsco',''),
                    'FDA Indication': ip.get('Indication_fda',''),
                })
            # Remove ID-matched entries from fuzzy matching pool
            matched_ids = set(id_pairs['RxCUI'])
            c = c[~c['RxCUI'].isin(matched_ids)]
            logging.info(f"Removed {len(matched_ids)} ID-matched CDSCO entries for fuzzy step")
        
        # Normalize drug names
        logging.info("Normalizing drug names")
        c["drug_norm"] = c["Drug Name"].apply(normalize)
        f["drug_norm"] = f["Drug Name"].apply(normalize)

        # Remove entries with empty normalized names
        initial_c, initial_f = len(c), len(f)
        c = c[c["drug_norm"] != ""]
        f = f[f["drug_norm"] != ""]
        logging.info(
            f"Dropped {initial_c - len(c)} CDSCO and {initial_f - len(f)} FDA entries with empty names"
        )

        # Drop duplicate normalized names
        c = c.drop_duplicates(subset=["drug_norm"])
        f = f.drop_duplicates(subset=["drug_norm"])
        logging.info(f"Unique normalized names: CDSCO={len(c)}, FDA={len(f)}")

        # Load synonyms mapping
        logging.info("Loading synonyms mapping (if any)")
        synonyms = load_synonyms()

        # Fuzzy matching with Jaccard pre-filtering
        logging.info("Starting fuzzy matching")
        for _, row in c.iterrows():
            # Get normalized name and apply synonyms mapping
            name = row["drug_norm"]
            name = synonyms.get(name, name)
            # Pre-filter candidates by Jaccard similarity
            candidates = f[
                f["drug_norm"].apply(lambda x: jaccard(name, x) >= jaccard_threshold)
            ]
            for _, fda_row in candidates.iterrows():
                # Compute similarity scores
                jw_score = jaro(name, fda_row["drug_norm"])
                token_score = fuzz.token_set_ratio(name, fda_row["drug_norm"])
                # Levenshtein ratio fallback
                ratio_score = fuzz.ratio(name, fda_row["drug_norm"])
                # Use consensus check for high-confidence matching
                if is_high_confidence_match(jw_score, token_score, ratio_score, threshold, token_threshold, ratio_threshold):
                    match_data = {
                        "CDSCO Drug Name": row["Drug Name"],
                        "FDA Drug Name": fda_row["Drug Name"],
                        "Similarity Score": jw_score,  # Jaro-Winkler
                        "Token Score": token_score,    # RapidFuzz token-set ratio
                        "Ratio Score": ratio_score,    # Levenshtein ratio
                        "CDSCO Approval Date": row.get("Date of Approval", ""),
                        "FDA Approval Date": fda_row.get("Date of Approval", ""),
                        "CDSCO Indication": row.get("Indication", ""),
                        "FDA Indication": fda_row.get("Indication", ""),
                    }
                    
                    # Add combination drug info if available
                    if use_exploded:
                        match_data["Original CDSCO Drug"] = row.get("Original Drug Name", row["Drug Name"])
                        match_data["Is Combination"] = row.get("Is Combination", False)
                    
                    matches.append(match_data)

        # Compile raw match results and emit only actual matches
        raw_df = pd.DataFrame(matches).drop_duplicates().reset_index(drop=True)
        # Ensure match type exists
        if 'Match Type' not in raw_df.columns:
            raw_df['Match Type'] = 'Fuzzy'
        else:
            raw_df['Match Type'] = raw_df['Match Type'].fillna('Fuzzy')
        # Select best match per CDSCO entry
        best_df = raw_df.sort_values(by='Similarity Score', ascending=False)
        best_df = best_df.drop_duplicates(subset=['CDSCO Drug Name'], keep='first').reset_index(drop=True)
        # Save only matched pairs
        best_df.to_csv(out_file, index=False)
        logging.info(f"Overlap results with only matched pairs written to {out_file}")
        logging.info(f"Found matches for {best_df['CDSCO Drug Name'].nunique()} CDSCO entries")
        logging.info(f"Total match records: {len(best_df)}")
    except Exception as e:
        logging.error(f"Error loading cleaned data: {e}")
        return


def main():
    """Parse command-line arguments and execute the drug matching process."""
    parser = argparse.ArgumentParser(
        description="Match CDSCO and FDA drug entries with high-confidence fuzzy matching"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=JARO_THRESHOLD,
        help="Jaro-Winkler similarity threshold",
    )
    parser.add_argument(
        "--jaccard-threshold",
        type=float,
        default=JACCARD_THRESHOLD,
        help="Jaccard similarity threshold for candidate filtering",
    )
    parser.add_argument(
        "--token-threshold",
        type=int,
        default=TOKEN_THRESHOLD,
        help="Token-set ratio threshold for fuzzy matching",
    )
    parser.add_argument(
        "--ratio-threshold",
        type=int,
        default=RATIO_THRESHOLD,
        help="Levenshtein ratio threshold for fuzzy matching",
    )
    parser.add_argument(
        "--out-file",
        type=str,
        default=str(PROC / "overlap.csv"),
        help="Output CSV file path",
    )
    parser.add_argument(
        "--use-exploded",
        action="store_true",
        help="Use the exploded CDSCO dataset with individual APIs",
    )
    args = parser.parse_args()
    run(
        threshold=args.threshold,
        jaccard_threshold=args.jaccard_threshold,
        token_threshold=args.token_threshold,
        ratio_threshold=args.ratio_threshold,
        out_file=args.out_file,
        use_exploded=args.use_exploded,
    )


if __name__ == "__main__":
    main()
