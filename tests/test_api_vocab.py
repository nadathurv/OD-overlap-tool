import pytest
import os
import pandas as pd
import pickle
from pathlib import Path
from src.utils.api_vocab import (
    extract_fda_vocabulary, 
    save_vocabulary, 
    load_vocabulary, 
    build_and_save_vocabulary
)
from src.utils.text import normalize
from src.config import PROC


@pytest.fixture
def mock_fda_dataframe():
    """Create a mock FDA dataframe for testing."""
    data = {
        "Drug Name": [
            "Paracetamol", 
            "Aspirin", 
            "Ibuprofen", 
            "Diclofenac Sodium",
            "Acetaminophen", # Same as Paracetamol but normalizes differently with current implementation
        ],
        "Indication": ["Pain", "Fever", "Inflammation", "Joint Pain", "Fever"]
    }
    df = pd.DataFrame(data)
    # Add normalized drug names
    df["drug_norm"] = df["Drug Name"].apply(normalize)
    return df


@pytest.fixture
def temp_vocab_path(tmp_path):
    """Create a temporary path for vocabulary storage."""
    return tmp_path / "test_vocab.pkl"


def test_extract_vocabulary(monkeypatch, mock_fda_dataframe):
    """Test extracting vocabulary from a DataFrame."""
    # Mock the read_parquet function to return our mock dataframe
    def mock_read_parquet(*args, **kwargs):
        return mock_fda_dataframe
    
    monkeypatch.setattr(pd, "read_parquet", mock_read_parquet)
    
    # Extract vocabulary
    vocab = extract_fda_vocabulary()
    
    # Check that we have the expected number of unique terms 
    # Note: With the current implementation, 'paracetamol' and 'acetaminophen' normalize to different terms
    assert len(vocab) == 5
    
    # Check that specific terms are in the vocabulary
    assert "paracetamol" in vocab
    assert "acetaminophen" in vocab
    assert "aspirin" in vocab
    assert "ibuprofen" in vocab
    assert "diclofenac sodium" in vocab


def test_save_and_load_vocabulary(temp_vocab_path):
    """Test saving and loading a vocabulary."""
    test_vocab = {"drug1", "drug2", "drug3"}
    
    # Save the vocabulary
    result = save_vocabulary(test_vocab, temp_vocab_path)
    assert result is True
    assert temp_vocab_path.exists()
    
    # Load the vocabulary
    loaded_vocab = load_vocabulary(temp_vocab_path)
    assert loaded_vocab == test_vocab


def test_load_nonexistent_vocabulary(temp_vocab_path):
    """Test loading a vocabulary that doesn't exist."""
    # Make sure the file doesn't exist
    if temp_vocab_path.exists():
        os.remove(temp_vocab_path)
    
    # Try to load a non-existent vocabulary
    vocab = load_vocabulary(temp_vocab_path)
    assert vocab == set() 