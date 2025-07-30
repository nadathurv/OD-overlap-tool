# Create a new test file for compare functions
import pytest
from src.analysis.compare import jaro


def test_jaro():
    assert jaro("hello", "hello") == 1.0
    assert jaro("hello", "hell") > 0.9
    assert jaro("hello", "world") < 0.5
    assert jaro("", "") == 1.0
    assert jaro("a", "") == 0.0
    assert jaro("", "a") == 0.0 