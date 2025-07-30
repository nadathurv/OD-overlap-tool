#!/usr/bin/env python3
"""
Unified drug overlap comparison pipeline.
Consolidates all comparison functionality into a single, efficient script.
"""

import pandas as pd
import argparse
from pathlib import Path
from typing import Dict, Optional, Tuple

from src.core.logging_config import setup_logging, get_logger
from src.core.matching import (
    optimized_fuzzy_matching, process_id_matches,
    jaro_winkler_similarity
)
from src.utils.text import normalize
from src.utils.synonyms import load_synonyms
from src.config import PROC

# Tightened thresholds based on empirical analysis
# Decision: Slightly increased Jaro threshold for higher precision
JARO_THRESHOLD = 0.72      # Tightened from 0.70 for better precision
JACCARD_THRESHOLD = 0.1    # Optimal for pre-filtering
TOKEN_THRESHOLD = 82       # Tightened from 80
RATIO_THRESHOLD = 82       # Tightened from 80

setup_logging()
logger = get_logger(__name__)


def load_dataset(name: str, raw_dir: Path) -> pd.DataFrame:
    """Load dataset variant with appropriate handling."""
    # Use processed parquet files when available
    if name in ["cdsco_clean", "fda_clean", "cdsco_exploded"]:
        path = PROC / f"{name}.parquet"
        if path.exists():
            return pd.read_parquet(path)
    
    # Load raw CSV variants
    csv_path = raw_dir / f"{name}.csv"
    
    if name == "cdsco_clean_names_only":
        # Names-only format: one drug name per line
        df = pd.read_csv(csv_path, header=None, names=['Drug Name'])
        df['Sr.No'] = range(1, len(df) + 1)
        df['Strength'] = df['Indication'] = df['Date of Approval'] = ''
    else:
        df = pd.read_csv(csv_path)
    
    return df


def enrich_names_only(names_df: pd.DataFrame, raw_dir: Path) -> pd.DataFrame:
    """
    Enrich names-only data with metadata from other sources.
    Decision: Use fuzzy matching with substring checks for robustness.
    """
    # Load enrichment sources
    unclean_df = load_dataset("unclean_cdsco", raw_dir)
    vikram_df = load_dataset("vikram_cdsco_clean", raw_dir)
    
    # Pre-compute normalized lookups for efficiency
    names_df['norm_name'] = names_df['Drug Name'].apply(normalize)
    
    # Build lookup dictionaries
    indication_lookup = {}
    for _, row in unclean_df.iterrows():
        norm = normalize(row.get('Drug Name', ''))
        if norm and pd.notna(row.get('Indication')):
            indication_lookup[norm] = row['Indication']
    
    date_lookup = {}
    for _, row in vikram_df.iterrows():
        norm = normalize(row.get('Drug Name', ''))
        if norm and pd.notna(row.get('Date of Approval')):
            date_lookup[norm] = row['Date of Approval']
    
    # Vectorized enrichment with fuzzy fallback
    def find_best_match(norm_name: str, lookup: Dict[str, str], 
                        threshold: float = 0.8) -> str:
        """Find best matching value using exact or fuzzy match."""
        # Exact match first
        if norm_name in lookup:
            return lookup[norm_name]
        
        # Fuzzy match with substring optimization
        best_score = 0
        best_value = ''
        
        for candidate_norm, value in lookup.items():
            # Quick substring check before expensive similarity
            if norm_name in candidate_norm or candidate_norm in norm_name:
                score = jaro_winkler_similarity(norm_name, candidate_norm)
                if score > best_score and score >= threshold:
                    best_score = score
                    best_value = value
        
        return best_value
    
    # Apply enrichment
    names_df['Indication'] = names_df['norm_name'].apply(
        lambda x: find_best_match(x, indication_lookup)
    )
    names_df['Date of Approval'] = names_df['norm_name'].apply(
        lambda x: find_best_match(x, date_lookup)
    )
    
    return names_df.drop(columns=['norm_name'])


def standardize_columns(df: pd.DataFrame, 
                       source: str = "fda") -> pd.DataFrame:
    """Standardize column names across different sources."""
    if source == "fda":
        # FDA uses different column names
        mappings = {
            'Generic Name': 'Drug Name',
            'Date Designated': 'Date of Approval',
            'Orphan Designation': 'Indication'
        }
        for old, new in mappings.items():
            if old in df.columns and new not in df.columns:
                df[new] = df[old]
    
    return df


def run_comparison(
    cdsco_source: str = "cdsco_clean",
    fda_source: str = "fda_clean",
    output_file: Optional[Path] = None,
    raw_dir: Path = Path("data/raw"),
    use_exploded: bool = False,
    **thresholds
) -> Tuple[pd.DataFrame, str]:
    """
    Core comparison logic with optimized performance.
    Returns (results_df, output_path) tuple.
    """
    # Use exploded variant if requested
    if use_exploded and cdsco_source == "cdsco_clean":
        cdsco_source = "cdsco_exploded"
    
    logger.info(f"Comparing {cdsco_source} vs {fda_source}")
    
    # Load and prepare datasets
    cdsco_df = load_dataset(cdsco_source, raw_dir)
    fda_df = load_dataset(fda_source, raw_dir)
    
    # Special handling for names-only variant
    if cdsco_source == "cdsco_clean_names_only":
        logger.info("Enriching names-only data")
        cdsco_df = enrich_names_only(cdsco_df, raw_dir)
    
    # Standardize columns
    fda_df = standardize_columns(fda_df, "fda")
    
    logger.info(
        f"Loaded {len(cdsco_df)} CDSCO and {len(fda_df)} FDA entries"
    )
    
    # Phase 1: ID-based matching (100% accurate)
    id_matches, remaining_cdsco = process_id_matches(cdsco_df, fda_df)
    logger.info(f"Found {len(id_matches)} RxNorm ID matches")
    
    # Phase 2: Prepare for fuzzy matching
    remaining_cdsco["drug_norm"] = remaining_cdsco["Drug Name"].apply(
        normalize
    )
    fda_df["drug_norm"] = fda_df["Drug Name"].apply(normalize)
    
    # Remove empty/invalid entries
    remaining_cdsco = remaining_cdsco[
        remaining_cdsco["drug_norm"].str.len() > 0
    ]
    fda_df = fda_df[fda_df["drug_norm"].str.len() > 0]
    
    # Deduplicate to ensure 1:1 matching
    remaining_cdsco = remaining_cdsco.drop_duplicates(subset=["drug_norm"])
    fda_df = fda_df.drop_duplicates(subset=["drug_norm"])
    
    logger.info(
        f"Fuzzy matching: {len(remaining_cdsco)} CDSCO vs {len(fda_df)} FDA"
    )
    
    # Load synonyms for name harmonization
    synonyms = load_synonyms()
    
    # Phase 3: Optimized fuzzy matching
    fuzzy_matches = optimized_fuzzy_matching(
        remaining_cdsco, fda_df, synonyms,
        jw_threshold=thresholds.get('jaro', JARO_THRESHOLD),
        jaccard_threshold=thresholds.get('jaccard', JACCARD_THRESHOLD),
        token_threshold=thresholds.get('token', TOKEN_THRESHOLD),
        ratio_threshold=thresholds.get('ratio', RATIO_THRESHOLD)
    )
    
    # Combine results
    all_matches = id_matches + fuzzy_matches
    results_df = pd.DataFrame(all_matches)
    
    # Post-processing: ensure 1:1 correspondence
    if not results_df.empty:
        # Remove duplicates, keeping highest scoring match
        results_df = results_df.sort_values(
            'Similarity Score', ascending=False
        )
        results_df = results_df.drop_duplicates(
            subset=['CDSCO Drug Name'], keep='first'
        )
        results_df = results_df.drop_duplicates(
            subset=['FDA Drug Name'], keep='first'
        )
        results_df = results_df.reset_index(drop=True)
    
    # Determine output path
    if output_file is None:
        suffix = f"_{cdsco_source}" if cdsco_source != "cdsco_clean" else ""
        output_file = PROC / f"overlap{suffix}.csv"
    
    # Save results
    output_file.parent.mkdir(parents=True, exist_ok=True)
    results_df.to_csv(output_file, index=False)
    
    logger.info(f"Saved {len(results_df)} matches to {output_file}")
    
    return results_df, str(output_file)


def main():
    parser = argparse.ArgumentParser(
        description="Unified drug overlap comparison pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Standard comparison
  python compare.py
  
  # Use exploded combinations
  python compare.py --use-exploded
  
  # Compare specific variant
  python compare.py --cdsco unclean_cdsco
  
  # Run all variants
  python compare.py --multi-variant
  
  # Custom thresholds
  python compare.py --jaro 0.75 --jaccard 0.15
        """
    )
    
    # Dataset selection
    parser.add_argument("--cdsco", default="cdsco_clean",
                        choices=["cdsco_clean", "cdsco_clean_names_only",
                                 "unclean_cdsco", "vikram_cdsco_clean"],
                        help="CDSCO dataset variant")
    parser.add_argument("--fda", default="fda_clean",
                        help="FDA dataset (default: fda_clean)")
    parser.add_argument("--use-exploded", action="store_true",
                        help="Use exploded combination drugs")
    
    # Output options
    parser.add_argument("--output", type=Path, help="Output CSV path")
    parser.add_argument("--raw-dir", type=Path, default=Path("data/raw"),
                        help="Raw data directory")
    
    # Multi-variant mode
    parser.add_argument("--multi-variant", action="store_true",
                        help="Run all CDSCO variants")
    
    # Threshold tuning
    parser.add_argument("--jaro", type=float, default=JARO_THRESHOLD,
                        help=f"Jaro-Winkler threshold (default: "
                             f"{JARO_THRESHOLD})")
    parser.add_argument("--jaccard", type=float, default=JACCARD_THRESHOLD,
                        help=f"Jaccard threshold (default: "
                             f"{JACCARD_THRESHOLD})")
    parser.add_argument("--token", type=int, default=TOKEN_THRESHOLD,
                        help=f"Token-set ratio threshold (default: "
                             f"{TOKEN_THRESHOLD})")
    parser.add_argument("--ratio", type=int, default=RATIO_THRESHOLD,
                        help=f"Levenshtein ratio threshold (default: "
                             f"{RATIO_THRESHOLD})")
    
    args = parser.parse_args()
    
    if args.multi_variant:
        # Run all CDSCO variants
        variants = ["cdsco_clean_names_only", "unclean_cdsco", 
                   "vikram_cdsco_clean"]
        print("\n🔄 Running multi-variant comparison...")
        
        for variant in variants:
            print(f"\n📊 Processing {variant}...")
            try:
                results_df, output_path = run_comparison(
                    cdsco_source=variant,
                    fda_source=args.fda,
                    raw_dir=args.raw_dir,
                    use_exploded=args.use_exploded,
                    jaro=args.jaro,
                    jaccard=args.jaccard,
                    token=args.token,
                    ratio=args.ratio
                )
                print(f"✅ {variant}: {len(results_df)} matches → "
                      f"{output_path}")
            except Exception as e:
                print(f"❌ {variant}: {e}")
        
        print("\n✨ Multi-variant comparison complete!")
    else:
        # Single comparison
        results_df, output_path = run_comparison(
            cdsco_source=args.cdsco,
            fda_source=args.fda,
            output_file=args.output,
            raw_dir=args.raw_dir,
            use_exploded=args.use_exploded,
            jaro=args.jaro,
            jaccard=args.jaccard,
            token=args.token,
            ratio=args.ratio
        )


if __name__ == "__main__":
    main() 