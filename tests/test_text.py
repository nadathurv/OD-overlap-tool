# Create a new test file for text utilities
import pytest
from src.utils.text import normalize, jaccard


def test_normalize():
    assert normalize("Hello, World!") == "hello world"
    assert normalize("  Multiple   spaces ") == "multiple spaces"
    assert normalize(None) == ""
    assert normalize("") == ""


def test_jaccard():
    assert jaccard("hello world", "hello") == 0.3333333333333333
    assert jaccard("hello world", "world hello") == 0.6666666666666666
    assert jaccard("", "") == 0.0
    assert jaccard("abc", "def") == 0.0
    assert jaccard("a b c", "a b c d") == 0.6
    assert jaccard("a b c", "") == 0.0 