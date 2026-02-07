
from dataclasses import dataclass

@dataclass
class CDSAnnotation:
    """Core CDS information and codon analysis"""
    transcript_id: str
    variant_id: str
    ref_cds_start: int
    ref_cds_stop: int
    ref_cds_seq: str
    ref_cds_len: int
    alt_cds_start: int
    alt_cds_stop: int
    alt_cds_seq: str
    alt_cds_len: int
    chromosome: str
    gene_id: str
    strand: str
    ref: str
    alt: str
    start_variant: int
    end_variant: int
    ref_cds_info: list
    alt_cds_info: list
    cds_in_transcript: bool
    ref_start_codon_pos: int = None
    ref_start_codon_exon: int = None
    ref_last_codon: str = None
    ref_valid_stop: bool = None
    ref_first_stop_codon: str = None
    ref_first_stop_pos: int = None
    ref_num_stop_codons: int = None
    ref_all_stop_codons: list = None
    ref_stop_codon_exons: list = None
    ref_is_premature: bool = None
    alt_start_codon_pos: int = None
    alt_start_codon_exon: int = None
    alt_last_codon: str = None
    alt_valid_stop: bool = None
    alt_first_stop_codon: str = None
    alt_first_stop_pos: int = None
    alt_num_stop_codons: int = None
    alt_all_stop_codons: list = None
    alt_stop_codon_exons: list = None
    alt_is_premature: bool = None
    start_loss: bool = None
    stop_loss: bool = None
    transcript_start: int = None
    transcript_end: int = None
    transcript_seq: str = None
    transcript_length: int = None
    alt_transcript_seq: str = None
    alt_transcript_length: int = None
    transcript_exon_info: list = None
    transcript_start_codon_pos: int = None
    transcript_start_codon_exon: int = None
    transcript_last_codon: str = None
    transcript_valid_stop: bool = None
    transcript_first_stop_codon: str = None
    transcript_first_stop_pos: int = None
    transcript_num_stop_codons: int = None
    transcript_all_stop_codons: list = None
    transcript_stop_codon_exons: list = None

@dataclass
class TranscriptFeatures:
    """Transcript-level features for NMD analysis"""
    utr3_length: int = None
    utr5_length: int = None
    total_exon_count: int = None
    upstream_exon_count: int = None
    downstream_exon_count: int = None
    ptc_to_start_codon: int = None
    ptc_less_than_150nt_to_start: bool = None
    ptc_exon_length: int = None
    stop_codon_distance: int = None
    ptc_to_intron: int = None
    likely_misannotated: bool = None

@dataclass
class NMDPrediction:
    """NMD escape rule evaluation"""
    nmd_last_exon_rule: bool = None
    nmd_50nt_penultimate_rule: bool = None
    nmd_long_exon_rule: bool = None
    nmd_start_proximal_rule: bool = None
    nmd_single_exon_rule: bool = None
    nmd_escape: bool = None
