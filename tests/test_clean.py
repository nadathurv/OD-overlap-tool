# Create a new test file for cleaning functions
import pytest
import pandas as pd
from src.data.clean import _norm_date, _norm_text, _strip_forms


def test_norm_date():
    assert _norm_date("2023") == "01/01/2023"
    assert _norm_date("12/2023") == "12/01/2023"
    assert _norm_date("01-02-23") == "02/01/2023"
    assert _norm_date(None) is pd.NaT


def test_norm_text():
    assert _norm_text("  Hello   World  ") == "Hello World"
    assert _norm_text(None) == ""


def test_strip_forms():
    assert _strip_forms("Paracetamol 500mg Tablet") == "Paracetamol"
    assert _strip_forms("Aspirin 100mg") == "Aspirin"
    assert _strip_forms("Cough Syrup") == "Cough Syrup" 