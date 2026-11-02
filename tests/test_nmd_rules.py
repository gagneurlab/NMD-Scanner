import pandas as pd
import pytest

from nmd_scanner.nmd_rules import nmd_efficiency_prediction_status
from nmd_scanner.schema import MODEL_INPUTS, PREDICTION_STATUSES


def _status_pair(**values):
    """A PTC variant that the model can predict: a new PTC and all model inputs set. ``values`` overrides columns."""

    pair = {name: 1 for name in MODEL_INPUTS}
    pair.update(
        unknown_reason=None,
        alt_has_ptc=True,
        ref_has_ptc=False,
        has_stop_codon=True,
        has_start_codon=True,
        start_loss=False,
    )
    pair.update(values)
    return pair


# The null model inputs of a variant-transcript pair with an unknown alt transcript: 13 of the 17 model inputs
_UNKNOWN_INPUTS = {
    name: None
    for name in MODEL_INPUTS
    if name not in ("total_exon_count", "utr5_length", "utr3_length", "transcript_length")
}


@pytest.mark.parametrize(
    ("values", "status"),
    [
        ({}, "ok"),
        ({"unknown_reason": "splice_site_destroyed", "alt_has_ptc": None, **_UNKNOWN_INPUTS}, "unknown_effect"),
        ({"alt_has_ptc": False}, "no_ptc"),
        ({"alt_has_ptc": None}, "no_ptc"),
        ({"ref_has_ptc": True}, "ref_ptc"),
        ({"has_stop_codon": False, "annotated_stop_distance": None, "utr3_length": None}, "no_annotated_stop"),
        ({"has_start_codon": False, "ptc_to_start_codon": None}, "no_annotated_start"),
        ({"start_loss": True, "ptc_to_start_codon": None}, "start_loss"),
        ({"start_loss": True}, "start_loss"),
        ({"ptc_to_start_codon": None}, "missing_input"),
        ({"ptc_to_exon_end": None}, "missing_input"),
        ({"ref_has_ptc": None}, "ok"),
        # overlapping reasons: the first one in PREDICTION_STATUSES wins
        (
            {"unknown_reason": "exon_boundary_ambiguous", "alt_has_ptc": None, "has_stop_codon": False},
            "unknown_effect",
        ),
        ({"alt_has_ptc": False, "ref_has_ptc": True, "ptc_to_start_codon": None}, "no_ptc"),
        ({"ref_has_ptc": True, "has_stop_codon": False, "utr3_length": None}, "ref_ptc"),
        ({"has_stop_codon": False, "utr3_length": None, "ptc_to_start_codon": None}, "no_annotated_stop"),
        (
            {"has_stop_codon": False, "has_start_codon": False, "utr3_length": None, "ptc_to_start_codon": None},
            "no_annotated_stop",
        ),
        ({"has_start_codon": False, "ptc_to_start_codon": None, "ptc_to_exon_end": None}, "no_annotated_start"),
        ({"start_loss": True, "ptc_to_start_codon": None, "ptc_to_exon_end": None}, "start_loss"),
    ],
)
def test_nmd_efficiency_prediction_status(values, status):
    table = pd.DataFrame([_status_pair(), _status_pair(**values)], index=[7, 3])

    result = nmd_efficiency_prediction_status(table)

    assert result.tolist() == ["ok", status]
    assert result.index.tolist() == [7, 3]
    assert result.dtype == pd.CategoricalDtype(PREDICTION_STATUSES)


def test_nmd_efficiency_prediction_status_values_are_the_documented_ones():
    assert PREDICTION_STATUSES == (
        "unknown_effect",
        "no_ptc",
        "ref_ptc",
        "no_annotated_stop",
        "no_annotated_start",
        "start_loss",
        "missing_input",
        "ok",
    )
