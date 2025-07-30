import pytest

from src.analysis.compare import is_high_confidence_match


def test_all_metrics_above_threshold():
    # All three metrics exceed thresholds
    assert is_high_confidence_match(0.9, 90, 90, 0.85, 85, 85)


def test_two_metrics_above_threshold():
    # Exactly two metrics exceed thresholds
    assert is_high_confidence_match(0.9, 80, 90, 0.85, 85, 85)  # jw and ratio
    assert is_high_confidence_match(0.8, 90, 90, 0.85, 85, 85)  # token and ratio
    assert is_high_confidence_match(0.9, 90, 80, 0.85, 85, 85)  # jw and token


def test_one_metric_above_threshold():
    # Only one metric exceeds threshold should fail
    assert not is_high_confidence_match(0.9, 80, 80, 0.85, 85, 85)
    assert not is_high_confidence_match(0.8, 90, 80, 0.85, 85, 85)
    assert not is_high_confidence_match(0.8, 80, 90, 0.85, 85, 85)


def test_edge_threshold_values():
    # Metrics equal to thresholds count as exceeding
    assert is_high_confidence_match(0.85, 85, 84, 0.85, 85, 85)
    assert is_high_confidence_match(0.84, 85, 85, 0.85, 85, 85) 