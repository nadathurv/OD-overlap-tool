#!/usr/bin/env python3
"""
synonyms_gen.py

Script to generate a synonyms.json file by fuzzy-matching unmatched CDSCO drug names to FDA canonical names.
"""
import pandas as pd
import json
from rapidfuzz import fuzz
from src.utils.text import normalize

# Load FDA cleaned dataset
fda = pd.read_parquet('data/processed/fda_clean.parquet')
# Build unique set of normalized FDA names
fda_norms = list({normalize(n) for n in fda['Drug Name'].dropna()})

# Load unmatched CDSCO names
with open('data/processed/unmatched_cdsco.txt', 'r', encoding='utf-8') as f:
    unmatched = [line.strip() for line in f if line.strip()]

# Generate synonym mappings
synonyms = {}
for name in unmatched:
    norm_name = normalize(name)
    best_score = 0
    best_match = None
    for fda_norm in fda_norms:
        score = fuzz.token_set_ratio(norm_name, fda_norm)
        if score > best_score:
            best_score = score
            best_match = fda_norm
    # Only keep confident mappings
    if best_score >= 85:
        synonyms[norm_name] = best_match

# Write to synonyms.json
with open('data/processed/synonyms.json', 'w', encoding='utf-8') as out:
    json.dump(synonyms, out, indent=2)

print(f"Wrote {len(synonyms)} synonym mappings to data/processed/synonyms.json") 