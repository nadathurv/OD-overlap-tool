#!/usr/bin/env python3
"""
synonyms_check.py

Script to verify which CDSCO raw names mapped via synonyms appear in overlap.csv.
"""
import json
import pandas as pd
from src.utils.text import normalize

# Load synonyms mapping (normalized variant -> normalized canonical)
with open('data/processed/synonyms.json', 'r') as f:
    synonyms = json.load(f)
# Load raw unmatched names
with open('data/processed/unmatched_cdsco.txt', 'r') as f:
    raw_names = [line.strip() for line in f if line.strip()]
# Build reverse map: norm -> list of raw names
rev_map = {}
for raw in raw_names:
    norm = normalize(raw)
    rev_map.setdefault(norm, []).append(raw)
# Load overlap results
ov = pd.read_csv('data/processed/overlap.csv')
matched = set(ov['CDSCO Drug Name'].unique())
# Check mapping results
matched_raw = []
unmatched_raw = []
for norm, canonical in synonyms.items():
    for raw in rev_map.get(norm, []):
        if raw in matched:
            matched_raw.append(raw)
        else:
            unmatched_raw.append(raw)
# Summarize
print(f"Total synonyms entries: {len(synonyms)}")
print(f"Raw names matched after synonyms: {len(matched_raw)}")
print(f"Raw names still unmatched: {len(unmatched_raw)}")

# Optionally list them
print("\nMatched via synonyms:")
for r in matched_raw:
    print(r)
print("\nUnmatched via synonyms:")
for r in unmatched_raw:
    print(r) 