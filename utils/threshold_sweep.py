#!/usr/bin/env python3
"""
threshold_sweep.py

Utility script to evaluate how many overlap matches we get under different Jaro-Winkler
and Jaccard threshold settings. This helps pick thresholds that maximize coverage.
"""
import pandas as pd
from itertools import product
from src.analysis.compare import run

# Define threshold ranges to sweep
JW_VALUES = [0.70, 0.80, 0.85, 0.90, 0.95]
JACCARD_VALUES = [0.10, 0.30, 0.50]

# Directory to write intermediate overlap files
OUT_DIR = "data/processed"


def sweep_thresholds(save_intermediate: bool = True):
    """
    Run the fuzzy matching pipeline across combinations of JW and Jaccard thresholds.
    Prints a summary table with total match counts and unique CDSCO matches.
    """
    results = []
    for jw, jacc in product(JW_VALUES, JACCARD_VALUES):
        out_file = f"{OUT_DIR}/overlap_jw{jw:.2f}_jac{jacc:.2f}.csv"
        print(f"Running thresholds: Jaro-W={jw:.2f}, Jaccard={jacc:.2f}")
        run(threshold=jw, jaccard_threshold=jacc, out_file=out_file)
        df = pd.read_csv(out_file)
        total_matches = len(df)
        cdsco_unique = df["CDSCO Drug Name"].nunique()
        fda_unique = df["FDA Drug Name"].nunique()
        results.append({
            "Jaro-Winkler": jw,
            "Jaccard": jacc,
            "TotalMatches": total_matches,
            "UniqueCDSCO": cdsco_unique,
            "UniqueFDA": fda_unique,
        })
    summary = pd.DataFrame(results)
    print("\nThreshold sweep summary:")
    print(summary.to_string(index=False))
    return summary


if __name__ == "__main__":
    sweep_thresholds() 