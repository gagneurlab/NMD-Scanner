"""
Tests for the new modular annotation API
"""
import pytest
from nmd_scanner.annotation import (
    annotate_cds,
    calculate_transcript_features,
    evaluate_nmd_escape,
    CDSAnnotation,
    TranscriptFeatures,
    NMDPrediction
)


def test_cds_annotation_model():
    """Test CDSAnnotation dataclass creation"""
    cds = CDSAnnotation(
        transcript_id="ENST00000001",
        variant_id="chr1:100:A>G",
        ref_cds_start=100,
        ref_cds_stop=200,
        ref_cds_seq="ATGCCC",
        ref_cds_len=6,
        alt_cds_start=100,
        alt_cds_stop=200,
        alt_cds_seq="ATGCGC",
        alt_cds_len=6,
        chromosome="chr1",
        gene_id="ENSG00000001",
        strand="+",
        ref="A",
        alt="G",
        start_variant=150,
        end_variant=151,
        ref_cds_info=[(1, 100)],
        alt_cds_info=[(1, 100)],
        cds_in_transcript=True,
        alt_is_premature=False
    )
    
    assert cds.transcript_id == "ENST00000001"
    assert cds.gene_id == "ENSG00000001"
    assert cds.variant_id == "chr1:100:A>G"
    assert cds.chromosome == "chr1"


def test_transcript_features_model():
    """Test TranscriptFeatures dataclass"""
    features = TranscriptFeatures(
        utr5_length=200,
        utr3_length=150,
        total_exon_count=5,
        upstream_exon_count=2,
        downstream_exon_count=2,
        ptc_exon_length=100
    )
    
    assert features.utr5_length == 200
    assert features.total_exon_count == 5


def test_nmd_prediction_model():
    """Test NMDPrediction dataclass"""
    prediction = NMDPrediction(
        nmd_last_exon_rule=False,
        nmd_50nt_penultimate_rule=False,
        nmd_long_exon_rule=False,
        nmd_start_proximal_rule=False,
        nmd_single_exon_rule=False,
        nmd_escape=False
    )
    
    assert prediction.nmd_escape == False
    assert prediction.nmd_last_exon_rule == False


def test_evaluate_nmd_escape_non_premature():
    """Test that non-premature variants don't trigger NMD"""
    cds = CDSAnnotation(
        transcript_id="ENST00000001",
        variant_id="chr1:100:A>G",
        ref_cds_start=100,
        ref_cds_stop=200,
        ref_cds_seq="ATGCCC",
        ref_cds_len=6,
        alt_cds_start=100,
        alt_cds_stop=200,
        alt_cds_seq="ATGCGC",
        alt_cds_len=6,
        chromosome="chr1",
        gene_id="ENSG00000001",
        strand="+",
        ref="A",
        alt="G",
        start_variant=150,
        end_variant=151,
        ref_cds_info=[(1, 100)],
        alt_cds_info=[(1, 100)],
        cds_in_transcript=True,
        alt_is_premature=False,
        alt_first_stop_pos=None,
        alt_start_codon_pos=0,
        transcript_exon_info=[(1, 300)]
    )
    
    features = TranscriptFeatures(
        total_exon_count=3,
        downstream_exon_count=1,
        ptc_exon_length=100
    )
    
    prediction = evaluate_nmd_escape(cds, features)
    
    assert prediction.nmd_escape == False
    assert prediction.nmd_last_exon_rule == False


def test_evaluate_nmd_escape_last_exon_rule():
    """Test NMD escape due to last exon rule"""
    cds = CDSAnnotation(
        transcript_id="ENST00000001",
        variant_id="chr1:100:A>G",
        ref_cds_start=100,
        ref_cds_stop=600,
        ref_cds_seq="ATG" * 100,
        ref_cds_len=300,
        alt_cds_start=100,
        alt_cds_stop=400,
        alt_cds_seq="ATG" * 50,
        alt_cds_len=150,
        chromosome="chr1",
        gene_id="ENSG00000001",
        strand="+",
        ref="A",
        alt="G",
        start_variant=250,
        end_variant=251,
        ref_cds_info=[(1, 100), (2, 100), (3, 100)],
        alt_cds_info=[(1, 100), (2, 50)],
        cds_in_transcript=True,
        alt_is_premature=True,
        alt_first_stop_pos=150,
        alt_start_codon_pos=0,
        alt_stop_codon_exons=[2],
        transcript_exon_info=[(1, 100), (2, 100), (3, 100)]
    )
    
    features = TranscriptFeatures(
        total_exon_count=3,
        downstream_exon_count=0,  # PTC in last exon
        ptc_exon_length=100
    )
    
    prediction = evaluate_nmd_escape(cds, features)
    
    assert prediction.nmd_last_exon_rule == True
    assert prediction.nmd_escape == True


def test_evaluate_nmd_escape_start_proximal_rule():
    """Test NMD escape due to start proximal rule"""
    cds = CDSAnnotation(
        transcript_id="ENST00000001",
        variant_id="chr1:100:A>G",
        ref_cds_start=100,
        ref_cds_stop=600,
        ref_cds_seq="ATG" * 100,
        ref_cds_len=300,
        alt_cds_start=100,
        alt_cds_stop=250,
        alt_cds_seq="ATG" * 50,
        alt_cds_len=150,
        chromosome="chr1",
        gene_id="ENSG00000001",
        strand="+",
        ref="A",
        alt="G",
        start_variant=150,
        end_variant=151,
        ref_cds_info=[(1, 300)],
        alt_cds_info=[(1, 150)],
        cds_in_transcript=True,
        alt_is_premature=True,
        alt_first_stop_pos=100,  # Within 150nt of start
        alt_start_codon_pos=0,
        alt_stop_codon_exons=[1],
        transcript_exon_info=[(1, 300)]
    )
    
    features = TranscriptFeatures(
        total_exon_count=1,
        downstream_exon_count=0,
        ptc_exon_length=300
    )
    
    prediction = evaluate_nmd_escape(cds, features)
    
    assert prediction.nmd_start_proximal_rule == True
    assert prediction.nmd_escape == True


def test_evaluate_nmd_escape_single_exon_rule():
    """Test NMD escape due to single exon rule"""
    cds = CDSAnnotation(
        transcript_id="ENST00000001",
        variant_id="chr1:100:A>G",
        ref_cds_start=100,
        ref_cds_stop=600,
        ref_cds_seq="ATG" * 100,
        ref_cds_len=300,
        alt_cds_start=100,
        alt_cds_stop=400,
        alt_cds_seq="ATG" * 100,
        alt_cds_len=300,
        chromosome="chr1",
        gene_id="ENSG00000001",
        strand="+",
        ref="A",
        alt="G",
        start_variant=250,
        end_variant=251,
        ref_cds_info=[(1, 300)],
        alt_cds_info=[(1, 300)],
        cds_in_transcript=True,
        alt_is_premature=True,
        alt_first_stop_pos=250,
        alt_start_codon_pos=0,
        alt_stop_codon_exons=[1],
        transcript_exon_info=[(1, 600)]
    )
    
    features = TranscriptFeatures(
        total_exon_count=1,  # Single exon transcript
        downstream_exon_count=0,
        ptc_exon_length=600
    )
    
    prediction = evaluate_nmd_escape(cds, features)
    
    assert prediction.nmd_single_exon_rule == True
    assert prediction.nmd_escape == True


def test_evaluate_nmd_escape_long_exon_rule():
    """Test NMD escape due to long exon rule (>407 nt)"""
    cds = CDSAnnotation(
        transcript_id="ENST00000001",
        variant_id="chr1:100:A>G",
        ref_cds_start=100,
        ref_cds_stop=1000,
        ref_cds_seq="ATG" * 300,
        ref_cds_len=900,
        alt_cds_start=100,
        alt_cds_stop=700,
        alt_cds_seq="ATG" * 200,
        alt_cds_len=600,
        chromosome="chr1",
        gene_id="ENSG00000001",
        strand="+",
        ref="A",
        alt="G",
        start_variant=400,
        end_variant=401,
        ref_cds_info=[(1, 300), (2, 600)],
        alt_cds_info=[(1, 300), (2, 300)],
        cds_in_transcript=True,
        alt_is_premature=True,
        alt_first_stop_pos=600,
        alt_start_codon_pos=0,
        alt_stop_codon_exons=[2],
        transcript_exon_info=[(1, 300), (2, 700)]  # Exon 2 is 700 nt
    )
    
    features = TranscriptFeatures(
        total_exon_count=2,
        downstream_exon_count=0,
        ptc_exon_length=700  # >407 nt
    )
    
    prediction = evaluate_nmd_escape(cds, features)
    
    assert prediction.nmd_long_exon_rule == True
    assert prediction.nmd_escape == True


def test_evaluate_nmd_escape_no_escape():
    """Test variant that triggers NMD (no escape)"""
    cds = CDSAnnotation(
        transcript_id="ENST00000001",
        variant_id="chr1:100:A>G",
        ref_cds_start=100,
        ref_cds_stop=600,
        ref_cds_seq="ATG" * 100,
        ref_cds_len=300,
        alt_cds_start=100,
        alt_cds_stop=400,
        alt_cds_seq="ATG" * 100,
        alt_cds_len=300,
        chromosome="chr1",
        gene_id="ENSG00000001",
        strand="+",
        ref="A",
        alt="G",
        start_variant=250,
        end_variant=251,
        ref_cds_info=[(1, 100), (2, 100), (3, 100)],
        alt_cds_info=[(1, 100), (2, 50)],
        cds_in_transcript=True,
        alt_is_premature=True,
        alt_first_stop_pos=200,  # >150 nt from start
        alt_start_codon_pos=0,
        alt_stop_codon_exons=[2],
        transcript_exon_info=[(1, 100), (2, 100), (3, 100), (4, 100)]
    )
    
    features = TranscriptFeatures(
        total_exon_count=4,
        downstream_exon_count=2,  # Not in last exon
        ptc_exon_length=100  # <407 nt
    )
    
    prediction = evaluate_nmd_escape(cds, features)
    
    assert prediction.nmd_last_exon_rule == False
    assert prediction.nmd_start_proximal_rule == False
    assert prediction.nmd_long_exon_rule == False
    assert prediction.nmd_single_exon_rule == False
    assert prediction.nmd_escape == False
