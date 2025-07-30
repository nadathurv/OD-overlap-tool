"""
Utilities for managing API vocabularies to support combination drug canonicalization.
"""

import json
import logging
import pandas as pd
from pathlib import Path
from typing import Set, List, Tuple
import pickle

from ..config import PROC, FDA_CLEAN

# Path for storing the vocabulary
FDA_VOCAB_PATH = PROC / "fda_api_vocab.pkl"

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")


def extract_fda_vocabulary() -> Set[str]:
    """
    Extract a vocabulary of normalized drug names from the cleaned FDA dataset.
    
    Returns:
        Set[str]: Set of normalized drug names that can be used as a reference vocabulary
    """
    try:
        logging.info(f"Loading FDA dataset from {FDA_CLEAN}")
        fda_df = pd.read_parquet(FDA_CLEAN)
        
        # Ensure drug_norm column exists
        if "drug_norm" not in fda_df.columns:
            from ..utils.text import normalize
            logging.info("Creating drug_norm column")
            fda_df["drug_norm"] = fda_df["Drug Name"].apply(normalize)
        
        # Extract unique normalized drug names
        vocab = set(fda_df["drug_norm"].dropna().unique())
        logging.info(f"Extracted {len(vocab)} unique API terms from FDA dataset")
        return vocab
        
    except Exception as e:
        logging.error(f"Error extracting FDA vocabulary: {e}")
        return set()


def save_vocabulary(vocab: Set[str], path: Path = FDA_VOCAB_PATH) -> bool:
    """
    Save the vocabulary to disk for later use.
    
    Args:
        vocab: Set of normalized drug names
        path: Path to save the vocabulary
        
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        path.parent.mkdir(exist_ok=True, parents=True)
        with open(path, 'wb') as f:
            pickle.dump(vocab, f)
        logging.info(f"Saved {len(vocab)} API terms to {path}")
        return True
    except Exception as e:
        logging.error(f"Error saving vocabulary: {e}")
        return False


def load_vocabulary(path: Path = FDA_VOCAB_PATH) -> Set[str]:
    """
    Load the vocabulary from disk.
    
    Args:
        path: Path to the saved vocabulary
        
    Returns:
        Set[str]: The loaded vocabulary or empty set if file not found
    """
    try:
        if not path.exists():
            logging.warning(f"Vocabulary file {path} not found")
            return set()
            
        with open(path, 'rb') as f:
            vocab = pickle.load(f)
        logging.info(f"Loaded {len(vocab)} API terms from {path}")
        return vocab
    except Exception as e:
        logging.error(f"Error loading vocabulary: {e}")
        return set()


def build_and_save_vocabulary() -> Set[str]:
    """
    Extract and save the FDA vocabulary in one step.
    
    Returns:
        Set[str]: The extracted vocabulary
    """
    vocab = extract_fda_vocabulary()
    if vocab:
        save_vocabulary(vocab)
    return vocab 