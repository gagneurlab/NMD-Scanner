# Import dependencies
import numpy as np
import pytest
from nmd_scanner.annotation.regressors import nmdeff_predict

def test_nmdeff_predict():
    """Test the NMD efficiency prediction model with sample data."""
    # Create sample input data (3 samples)
    n_samples = 3
    
    # Binary features
    start_loss = np.array([0, 1, 0])
    stop_loss = np.array([1, 0, 0])
    ptc_less_than_150nt_to_start = np.array([0, 1, 0])
    nmd_long_exon_rule = np.array([1, 0, 0])
    nmd_start_proximal_rule = np.array([0, 0, 1])
    nmd_single_exon_rule = np.array([0, 0, 0])
    nmd_escape = np.array([1, 1, 0])
    nmd_last_exon_rule = np.array([0, 1, 0])
    nmd_50nt_penultimate_rule = np.array([0, 0, 0])
    
    # Numeric features
    total_exon_count = np.array([10, 5, 8])
    downstream_exon_count = np.array([3, 0, 2])
    upstream_exon_count = np.array([7, 5, 6])
    ptc_to_start_codon = np.array([500, 200, 350])
    stop_codon_distance = np.array([150, 50, 100])
    ptc_exon_length = np.array([200, 150, 180])
    ptc_to_intron = np.array([80, 30, 60])
    utr5_length = np.array([100, 150, 120])
    utr3_length = np.array([300, 250, 280])
    transcript_length = np.array([2000, 1500, 1800])

    # Get predictions
    predictions = nmdeff_predict(
        start_loss, stop_loss, total_exon_count, ptc_less_than_150nt_to_start,
        nmd_long_exon_rule, nmd_start_proximal_rule, nmd_single_exon_rule, nmd_escape,
        downstream_exon_count, nmd_last_exon_rule, ptc_to_start_codon,
        stop_codon_distance, ptc_exon_length, ptc_to_intron, upstream_exon_count,
        nmd_50nt_penultimate_rule, utr5_length, utr3_length, transcript_length
    )
    
    # Assertions
    assert predictions is not None, "Predictions should not be None"
    assert len(predictions) == n_samples, f"Expected {n_samples} predictions, got {len(predictions)}"
    assert isinstance(predictions, np.ndarray), "Predictions should be a numpy array"
    assert predictions.dtype in [np.float32, np.float64], "Predictions should be float values"
    
    # Check that predictions are numeric (can be any real number for regression)
    assert np.all(np.isfinite(predictions)), "All predictions should be finite numbers"


def test_nmdeff_predict_single_sample():
    """Test the NMD efficiency prediction with a single sample."""
    # Single sample input
    predictions = nmdeff_predict(
        np.array([0]), np.array([0]), np.array([10]), np.array([0]),
        np.array([0]), np.array([0]), np.array([0]), np.array([1]),
        np.array([3]), np.array([0]), np.array([500]), np.array([150]),
        np.array([200]), np.array([80]), np.array([7]), np.array([0]),
        np.array([100]), np.array([300]), np.array([2000])
    )
    
    assert len(predictions) == 1, "Should return exactly one prediction"
    assert isinstance(predictions[0], (np.floating, float)), "Prediction should be a float"
    assert np.isfinite(predictions[0]), "Prediction should be a finite number"

