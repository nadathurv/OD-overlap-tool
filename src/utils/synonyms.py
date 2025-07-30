"""
Utilities for loading drug name synonyms to improve matching accuracy.
"""
import json
import logging
from pathlib import Path
from ..config import ROOT

# Path to a JSON file containing synonym mappings: {"variant_norm": "canonical_norm"}
SYNONYMS_PATH = ROOT / "data" / "processed" / "synonyms.json"

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

def load_synonyms():
    """
    Load synonyms from SYNONYMS_PATH.

    Returns:
        dict: mapping normalized variant -> normalized canonical name
    """
    if not SYNONYMS_PATH.exists():
        logging.warning(f"Synonyms file not found at {SYNONYMS_PATH}; proceeding without synonyms.")
        return {}
    try:
        with open(SYNONYMS_PATH, 'r', encoding='utf-8') as f:
            data = json.load(f)
        # Ensure keys and values are lowercase normalized
        norm_map = {k.lower(): v.lower() for k, v in data.items()}
        logging.info(f"Loaded {len(norm_map)} synonyms from {SYNONYMS_PATH}")
        return norm_map
    except Exception as e:
        logging.error(f"Error loading synonyms: {e}")
        return {} 