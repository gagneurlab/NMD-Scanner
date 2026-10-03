# Import dependencies
import gzip
import re
from pathlib import Path

import pandas as pd
import pyranges as pr
import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

import nmd_scanner
from nmd_scanner.cli import main
from nmd_scanner.scan import merge_stop_codons_into_cds

# pytest-fixtures as inputs for the tests


@pytest.fixture(scope="session")
def gtf_path():
    return "resources/chr18.gtf.gz"


@pytest.fixture(scope="session")
def vcf_path():
    return "resources/part-00241-61a0abbf-fbf9-444f-8287-4e46ad4b9b7b-c000.vcf"


@pytest.fixture(scope="session")
def fasta_path():
    return "resources/chr18.fa.gz"


# Create the test functions


# Test reading VCF file
def test_read_vcf_file(vcf_path):

    gr = nmd_scanner.scan.read_vcf(vcf_path)
    assert gr is not None
    assert gr.df.shape[0] > 0
    assert "Chromosome" in gr.df.columns
    assert "Start" in gr.df.columns
    assert "End" in gr.df.columns

    print(gr.df.head())
    print(gr.df.shape)


def test_read_vcf_rejects_multiallelic(tmp_path):
    vcf = tmp_path / "multiallelic.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\tv1\tA\tT\t.\t.\t.\n"
        "chr1\t200\tv2\tC\tG,GT\t.\t.\t.\n"
    )
    with pytest.raises(ValueError, match="multi-allelic"):
        nmd_scanner.scan.read_vcf(str(vcf))


def test_read_vcf_accepts_single_allelic(tmp_path):
    vcf = tmp_path / "single.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\tv1\tA\tT\t.\t.\t.\n"
        "chr1\t200\tv2\tC\tG\t.\t.\t.\n"
    )
    gr = nmd_scanner.scan.read_vcf(str(vcf))
    assert gr.df.shape[0] == 2


def test_read_vcf_keeps_text_fields_as_written(tmp_path):
    vcf = tmp_path / "text_fields.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "22\t100\t12345\tA\tNA\t.\tPASS\tNA\n"
        "22\t200\tNA\tC\tG\t50\t.\t.\n"
    )
    df = nmd_scanner.scan.read_vcf(str(vcf)).df
    assert (df["Chromosome"] == "22").all()
    assert df["ID"].tolist() == ["12345", "NA"]
    assert df["Alt"].tolist() == ["NA", "G"]
    assert df["Qual"].tolist() == [".", "50"]
    assert df["Filter"].tolist() == ["PASS", "."]
    assert df["Info"].tolist() == ["NA", "."]
    assert df["Start"].tolist() == [99, 199]
    assert df["End"].tolist() == [100, 200]


# Test reading GTF file
def test_read_gtf_file(gtf_path):
    gr = nmd_scanner.scan.read_gtf(gtf_path)
    assert gr is not None
    assert gr.df.shape[0] > 0
    assert "Chromosome" in gr.df.columns
    print(gr.df.head())


def _coding_rows(rows):
    """CDS and stop_codon rows of one transcript: (Feature, exon_number, Start, End, Strand)."""
    return pd.DataFrame(
        [
            {"transcript_id": "tx", "Feature": f, "exon_number": str(e), "Start": s, "End": en, "Strand": st}
            for f, e, s, en, st in rows
        ]
    )


def _intervals(df):
    return sorted(zip(df["exon_number"], df["Start"], df["End"], df["Feature"]))


def test_merge_stop_codons_into_cds_gtf():
    # GTF: the stop codon follows the last CDS on + strand, precedes it on - strand
    plus = _coding_rows([("CDS", 1, 100, 150, "+"), ("CDS", 2, 200, 250, "+"), ("stop_codon", 2, 250, 253, "+")])
    assert _intervals(merge_stop_codons_into_cds(plus)) == [(1, 100, 150, "CDS"), (2, 200, 253, "CDS")]

    minus = _coding_rows([("CDS", 2, 500, 550, "-"), ("CDS", 1, 800, 850, "-"), ("stop_codon", 2, 497, 500, "-")])
    assert _intervals(merge_stop_codons_into_cds(minus)) == [(1, 800, 850, "CDS"), (2, 497, 550, "CDS")]


def test_merge_stop_codons_into_cds_without_stop_codon():
    # cds_end_NF: no stop_codon row, so the CDS stays as it is
    rows = _coding_rows([("CDS", 1, 100, 150, "+"), ("CDS", 2, 200, 250, "+")])
    assert _intervals(merge_stop_codons_into_cds(rows)) == [(1, 100, 150, "CDS"), (2, 200, 250, "CDS")]


def test_merge_stop_codons_into_cds_split_stop_codon():
    # stop codon split across an intron: 2 bases at the end of exon 2, 1 base at the start of exon 3,
    # which has no CDS row
    rows = _coding_rows(
        [
            ("CDS", 1, 100, 150, "+"),
            ("CDS", 2, 200, 248, "+"),
            ("stop_codon", 2, 248, 250, "+"),
            ("stop_codon", 3, 300, 301, "+"),
        ]
    )
    assert _intervals(merge_stop_codons_into_cds(rows)) == [
        (1, 100, 150, "CDS"),
        (2, 200, 250, "CDS"),
        (3, 300, 301, "CDS"),
    ]


def test_merge_stop_codons_into_cds_split_stop_codon_minus_strand():
    # stop codon split across an intron on the - strand: 2 bases at the lower end of exon 2, 1 base in exon 3,
    # which lies at lower coordinates and has no CDS row
    rows = _coding_rows(
        [
            ("CDS", 1, 800, 850, "-"),
            ("CDS", 2, 500, 550, "-"),
            ("stop_codon", 2, 498, 500, "-"),
            ("stop_codon", 3, 400, 401, "-"),
        ]
    )
    assert _intervals(merge_stop_codons_into_cds(rows)) == [
        (1, 800, 850, "CDS"),
        (2, 498, 550, "CDS"),
        (3, 400, 401, "CDS"),
    ]


def test_merge_stop_codons_into_cds_has_stop_codon():
    # the flag is per transcript and comes from the stop_codon rows
    with_stop = _coding_rows([("CDS", 1, 100, 150, "+"), ("CDS", 2, 200, 250, "+"), ("stop_codon", 2, 250, 253, "+")])
    without_stop = _coding_rows([("CDS", 1, 100, 150, "+"), ("CDS", 2, 200, 250, "+")]).assign(transcript_id="tx_nf")
    merged = merge_stop_codons_into_cds(pd.concat([with_stop, without_stop], ignore_index=True))
    assert merged.groupby("transcript_id")["has_stop_codon"].agg(set).to_dict() == {"tx": {True}, "tx_nf": {False}}


def test_merge_stop_codons_into_cds_warns_without_stop_codon_rows(caplog):
    rows = _coding_rows([("CDS", 1, 100, 150, "+"), ("CDS", 2, 200, 250, "+")])
    with caplog.at_level("WARNING", logger="nmd_scanner.scan"):
        merge_stop_codons_into_cds(rows)
    assert "No stop_codon rows" in caplog.text
    assert "the merge needs the stop_codon rows of the GTF" in caplog.text

    caplog.clear()
    with_stop = _coding_rows([("CDS", 1, 100, 150, "+"), ("stop_codon", 1, 150, 153, "+")])
    with caplog.at_level("WARNING", logger="nmd_scanner.scan"):
        merge_stop_codons_into_cds(with_stop)
    assert "No stop_codon rows" not in caplog.text


def test_merge_stop_codons_into_cds_rejects_gap():
    rows = _coding_rows([("CDS", 1, 100, 150, "+"), ("stop_codon", 1, 160, 163, "+")])
    with pytest.raises(ValueError, match="tx"):
        merge_stop_codons_into_cds(rows)


def test_merge_stop_codons_into_cds_rejects_gap_minus_strand():
    # on the - strand the stop codon lies below the CDS: stop codon [496, 499) and CDS [500, 550) leave out base 499
    rows = _coding_rows([("CDS", 1, 500, 550, "-"), ("stop_codon", 1, 496, 499, "-")])
    with pytest.raises(ValueError, match="tx"):
        merge_stop_codons_into_cds(rows)


def _coding_sequence(coding, fasta):
    coding = coding.sort_values("Start")
    seq = "".join(fasta[c][s:e].seq.upper() for c, s, e in zip(coding["Chromosome"], coding["Start"], coding["End"]))
    return str(Seq(seq).reverse_complement()) if coding["Strand"].iloc[0] == "-" else seq


def test_merge_stop_codons_into_cds_on_real_transcripts():
    gtf_df = nmd_scanner.scan.read_gtf("resources/chr18.gtf.gz").df
    fasta = Fasta("resources/chr18.fa.gz")
    rows = gtf_df[gtf_df["Feature"].isin(["CDS", "stop_codon"])]

    # ENST00000399496.8: stop codon split across an intron (two stop_codon rows)
    split = rows[rows["transcript_id"] == "ENST00000399496.8"]
    assert (split["Feature"] == "stop_codon").sum() == 2
    seq = _coding_sequence(merge_stop_codons_into_cds(split), fasta)
    assert len(seq) % 3 == 0
    assert seq[-3:] in {"TAA", "TAG", "TGA"}

    # ENST00000454642.3: minus strand, stop codon split across an intron (two stop_codon rows)
    split_minus = rows[rows["transcript_id"] == "ENST00000454642.3"]
    assert (split_minus["Feature"] == "stop_codon").sum() == 2
    assert split_minus["Strand"].iloc[0] == "-"
    seq = _coding_sequence(merge_stop_codons_into_cds(split_minus), fasta)
    assert len(seq) % 3 == 0
    assert seq[-3:] in {"TAA", "TAG", "TGA"}

    # a cds_end_NF transcript has no stop codon: its CDS rows stay as they are
    cds_end_nf = gtf_df.loc[(gtf_df["Feature"] == "transcript") & gtf_df["tag"].str.contains("cds_end_NF", na=False)]
    tx = cds_end_nf["transcript_id"].iloc[0]
    cds = rows[rows["transcript_id"] == tx]
    assert (cds["Feature"] == "stop_codon").sum() == 0
    merged = merge_stop_codons_into_cds(cds)
    assert sorted(zip(merged["Start"], merged["End"])) == sorted(zip(cds["Start"], cds["End"]))


# Test reading FASTA file


def test_read_fasta_file(fasta_path):
    fasta = nmd_scanner.scan.read_fasta(fasta_path)
    assert fasta is not None

    keys = list(fasta.keys())
    assert isinstance(keys, list)
    assert len(keys) > 0


def test_compute_exon_numbers():

    # On + Strand: Smallest exon number is the Start, Largest exon number is the End.
    df1 = pd.DataFrame(
        [
            ["chr1", 100, 200, "+", "exon", "TX1", "G1"],  # exon → 1
            ["chr1", 300, 400, "+", "exon", "TX1", "G1"],  # exon → 2
            ["chr1", 320, 400, "+", "CDS", "TX1", "G1"],  # CDS on exon 2
        ],
        columns=["Chromosome", "Start", "End", "Strand", "Feature", "transcript_id", "gene_id"],
    )
    out1 = nmd_scanner.compute_exon_numbers(pr.PyRanges(df1)).df

    tx1_exons = out1[(out1.Feature == "exon") & (out1.transcript_id == "TX1")].sort_values("Start")
    assert list(tx1_exons["exon_number"]) == [1, 2]
    tx1_cds = out1[(out1.Feature == "CDS") & (out1.transcript_id == "TX1")].iloc[0]
    assert tx1_cds["exon_number"] == 2

    # On - Strand: Smallest exon number is the Start, Largest exon number is the End.
    df2 = pd.DataFrame(
        [
            ["chr1", 100, 200, "-", "exon", "TX2", "G2"],  # exon_number → 2 (reverse order)
            ["chr1", 300, 400, "-", "exon", "TX2", "G2"],  # exon_number → 1
            ["chr1", 120, 180, "-", "CDS", "TX2", "G2"],  # CDS on exon 2
        ],
        columns=["Chromosome", "Start", "End", "Strand", "Feature", "transcript_id", "gene_id"],
    )

    out2 = nmd_scanner.compute_exon_numbers(pr.PyRanges(df2)).df

    tx2_exons = out2[(out2.Feature == "exon") & (out2.transcript_id == "TX2")].sort_values("Start")
    assert list(tx2_exons["exon_number"]) == [2, 1]
    tx2_cds = out2[(out2.Feature == "CDS") & (out2.transcript_id == "TX2")].iloc[0]
    assert tx2_cds["exon_number"] == 2

    # multiple CDS sequences on minus strand
    df4 = pd.DataFrame(
        [
            # exons
            ["chr1", 30, 50, "-", "exon", "TX2b", "G2b"],  # exon number 4
            ["chr1", 100, 200, "-", "exon", "TX2b", "G2b"],  # exon number 3
            ["chr1", 300, 350, "-", "exon", "TX2b", "G2b"],  # exon number 2
            ["chr1", 500, 600, "-", "exon", "TX2b", "G2b"],  # exon number 1
            # cds segments
            ["chr1", 500, 590, "-", "CDS", "TX2b", "G2b"],  # exon_number 1
            ["chr1", 300, 350, "-", "CDS", "TX2b", "G2b"],  # exon_number 2
            ["chr1", 130, 200, "-", "CDS", "TX2b", "G2b"],  # exon_number 3
        ],
        columns=["Chromosome", "Start", "End", "Strand", "Feature", "transcript_id", "gene_id"],
    )
    out4 = nmd_scanner.compute_exon_numbers(pr.PyRanges(df4)).df

    # check exons exon-numbers
    tx2b_exons = out4[(out4.Feature == "exon") & (out4.transcript_id == "TX2b")].sort_values("Start")
    assert list(tx2b_exons["exon_number"]) == [4, 3, 2, 1]

    # check CDS exon-numbers (coordinates must match the input rows)
    assert (
        int(
            out4[(out4.Feature == "CDS") & (out4.transcript_id == "TX2b") & (out4.Start == 500) & (out4.End == 590)][
                "exon_number"
            ].iloc[0]
        )
        == 1
    )
    assert (
        int(
            out4[(out4.Feature == "CDS") & (out4.transcript_id == "TX2b") & (out4.Start == 300) & (out4.End == 350)][
                "exon_number"
            ].iloc[0]
        )
        == 2
    )
    assert (
        int(
            out4[(out4.Feature == "CDS") & (out4.transcript_id == "TX2b") & (out4.Start == 130) & (out4.End == 200)][
                "exon_number"
            ].iloc[0]
        )
        == 3
    )

    # two different transripts (should be numbered independently)
    df3 = pd.DataFrame(
        [
            # TXA (+)
            ["chr1", 100, 150, "+", "exon", "TXA", "GA"],  # → 1
            ["chr1", 200, 250, "+", "exon", "TXA", "GA"],  # → 2
            # TXB (+)
            ["chr1", 500, 600, "+", "exon", "TXB", "GB"],  # → 1
            ["chr1", 700, 800, "+", "exon", "TXB", "GB"],  # → 2
            ["chr1", 900, 1000, "+", "exon", "TXB", "GB"],  # → 3
        ],
        columns=["Chromosome", "Start", "End", "Strand", "Feature", "transcript_id", "gene_id"],
    )
    out3 = nmd_scanner.compute_exon_numbers(pr.PyRanges(df3)).df

    ex_txA = out3[(out3.Feature == "exon") & (out3.transcript_id == "TXA")].sort_values("Start")
    ex_txB = out3[(out3.Feature == "exon") & (out3.transcript_id == "TXB")].sort_values("Start")
    assert list(ex_txA["exon_number"]) == [1, 2]
    assert list(ex_txB["exon_number"]) == [1, 2, 3]

    # TODO: maybe add the edge case if:
    # 1. CDS does not overlap any exon --> should not crash but exon_number should stay missing
    # 2. CDS overlaps two exons --> should it inherit the exon_number with the maximum overlap??


def test_compute_exon_numbers_with_str_exon_number_column():
    """
    A GTF read from file has a ``str`` exon_number column (pandas 3), with missing values on
    features that have no exon number. Computed exon numbers must be ints, not written into
    the str column.
    """

    df = pd.DataFrame(
        [
            ["chr1", 0, 500, "+", "gene", "TX1", "G1", None],
            ["chr1", 100, 200, "+", "exon", "TX1", "G1", "7"],
            ["chr1", 300, 400, "+", "exon", "TX1", "G1", "8"],
            ["chr1", 320, 400, "+", "CDS", "TX1", "G1", "8"],
        ],
        columns=["Chromosome", "Start", "End", "Strand", "Feature", "transcript_id", "gene_id", "exon_number"],
    )
    df["exon_number"] = df["exon_number"].astype("str")
    out = nmd_scanner.compute_exon_numbers(pr.PyRanges(df)).df

    exons = out[out.Feature == "exon"].sort_values("Start")
    assert list(exons["exon_number"]) == [1, 2]
    cds = out[out.Feature == "CDS"].iloc[0]
    assert cds["exon_number"] == 2
    assert not isinstance(cds["exon_number"], str)
    assert pd.isna(out[out.Feature == "gene"].iloc[0]["exon_number"])
    assert pd.api.types.is_integer_dtype(out["exon_number"])


# Test annotation format detection and dispatch


def test_detect_annotation_format():
    assert nmd_scanner.scan.detect_annotation_format("annotation.gtf") == "gtf"
    assert nmd_scanner.scan.detect_annotation_format("annotation.gtf.gz") == "gtf"
    assert nmd_scanner.scan.detect_annotation_format("annotation.gff3") == "gff3"
    assert nmd_scanner.scan.detect_annotation_format("annotation.gff3.gz") == "gff3"
    assert nmd_scanner.scan.detect_annotation_format("annotation.gff") == "gff3"
    assert nmd_scanner.scan.detect_annotation_format("annotation.gff.gz") == "gff3"
    assert nmd_scanner.scan.detect_annotation_format("ANNOTATION.GTF") == "gtf"

    with pytest.raises(ValueError, match="Cannot detect annotation format"):
        nmd_scanner.scan.detect_annotation_format("annotation.txt")


def _sorted_rows(df):
    df = df.copy()
    for col in df.columns:
        if isinstance(df[col].dtype, pd.CategoricalDtype):
            df[col] = df[col].astype(str)
    return df.sort_values(["transcript_id", "Feature", "Start", "End"]).reset_index(drop=True)


def test_read_annotation_gives_the_coding_regions_of_a_gtf(gtf_path):
    """The CDS rows are the coding regions; the other rows, minus the stop_codon rows, are the GTF rows."""
    via_annotation = nmd_scanner.scan.read_annotation(gtf_path).df
    gtf_df = nmd_scanner.scan.read_gtf(gtf_path).df
    assert not (via_annotation["Feature"] == "stop_codon").any()

    is_cds = via_annotation["Feature"] == "CDS"
    coding = merge_stop_codons_into_cds(gtf_df).astype({"exon_number": "Int64", "has_stop_codon": "boolean"})
    pd.testing.assert_frame_equal(_sorted_rows(via_annotation[is_cds]), _sorted_rows(coding))

    other = gtf_df[~gtf_df["Feature"].isin(["CDS", "stop_codon"])].astype({"exon_number": "Int64"})
    expected_other = other.assign(has_stop_codon=pd.Series(pd.NA, index=other.index, dtype="boolean"))
    pd.testing.assert_frame_equal(_sorted_rows(via_annotation[~is_cds]), _sorted_rows(expected_other))


def test_read_annotation_gives_the_hand_derived_coding_regions_of_the_fixture_gtf(tmp_path):
    """CDS plus stop_codon bases per exon of _GENCODE_GTF, 0-based half-open."""
    df = nmd_scanner.scan.read_annotation(_write(tmp_path, "gencode.gtf", _GENCODE_GTF)).df
    cds = df[df["Feature"] == "CDS"]
    got = sorted(
        zip(cds["transcript_id"].astype(str), cds["Start"], cds["End"], cds["exon_number"], cds["has_stop_codon"])
    )
    assert got == [
        ("ENST001.1", 1050, 1200, 1, True),
        ("ENST001.1", 1499, 2000, 2, True),
        ("ENST002.1", 3999, 4300, 2, True),
        ("ENST002.1", 5099, 5299, 1, True),
        ("ENST003.1", 7050, 7100, 1, True),
        # second base of the split stop codon, in an exon without CDS
        ("ENST003.1", 7199, 7200, 2, True),
        ("ENST004.1", 8049, 8300, 1, False),
        ("ENST005.1", 9050, 9098, 1, False),
        ("ENST006.1", 100, 160, 1, False),
        ("ENST007.1", 9801, 9849, 1, False),
    ]


def test_read_annotation_reassigns_gtf_exon_numbers_before_the_merge(tmp_path):
    """A GTF without exon_number attribute works with reassign_exons, also with a split stop codon (ENST003.1)."""
    without_numbers = re.sub(r" exon_number \d+;", "", _GENCODE_GTF)
    assert "exon_number" not in without_numbers
    path = _write(tmp_path, "no_numbers.gtf", without_numbers)

    reassigned = nmd_scanner.scan.read_annotation(path, reassign_exons=True)
    as_given = nmd_scanner.scan.read_annotation(_write(tmp_path, "gencode.gtf", _GENCODE_GTF))
    pd.testing.assert_frame_equal(_cds_exon_table(reassigned), _cds_exon_table(as_given))


def _write(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(content)
    return str(path)


# A tiny GENCODE-style GTF/GFF3 pair, and the equivalent Ensembl-style GTF/GFF3 pair, for the
# same loci: a + and a - strand transcript, a stop codon split across an intron (ENST003.1), a
# transcript without stop codon (cds_end_NF, ENST004.1) whose last 3 bases read TAA out of frame,
# and a chrM transcript whose CDS ends in AGA (ENST006.1). The GFF3 CDS includes the stop codon,
# the GTF CDS does not; the coding regions (CDS plus stop codon) are the same. The bases are in
# _FIXTURE_BASES.
# Ensembl: the GFF3 has no stop_codon rows, read_gff3 takes has_stop_codon from the FASTA. The GTF
# has stop codon rows for AGA on MT. ENSTE004 has ensembl_end_phase 2, so it has no stop codon.
# GENCODE: the GTF has no stop codon for AGA on chrM. ENST005.1 (+ strand) and ENST007.1 (- strand,
# split across an intron) are tagged cds_end_NF, but their CDS ends in a complete stop codon without
# stop_codon rows; the GENCODE GTF has these 3 bases as UTR.

_GENCODE_GTF = """\
chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSG001.1"; gene_type "protein_coding"; gene_name "GENE1";
chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "ENSG001.1"; transcript_id "ENST001.1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tgene_id "ENSG001.1"; transcript_id "ENST001.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chr1\tHAVANA\texon\t1500\t2000\t.\t+\t.\tgene_id "ENSG001.1"; transcript_id "ENST001.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 2;
chr1\tHAVANA\tCDS\t1051\t1200\t.\t+\t0\tgene_id "ENSG001.1"; transcript_id "ENST001.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chr1\tHAVANA\tCDS\t1500\t1997\t.\t+\t0\tgene_id "ENSG001.1"; transcript_id "ENST001.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 2;
chr1\tHAVANA\tstop_codon\t1998\t2000\t.\t+\t0\tgene_id "ENSG001.1"; transcript_id "ENST001.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 2;
chr1\tHAVANA\tgene\t4000\t5300\t.\t-\t.\tgene_id "ENSG002.1"; gene_type "protein_coding"; gene_name "GENE2";
chr1\tHAVANA\ttranscript\t4000\t5300\t.\t-\t.\tgene_id "ENSG002.1"; transcript_id "ENST002.1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1\tHAVANA\texon\t5100\t5300\t.\t-\t.\tgene_id "ENSG002.1"; transcript_id "ENST002.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chr1\tHAVANA\texon\t4000\t4300\t.\t-\t.\tgene_id "ENSG002.1"; transcript_id "ENST002.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 2;
chr1\tHAVANA\tCDS\t5100\t5299\t.\t-\t0\tgene_id "ENSG002.1"; transcript_id "ENST002.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chr1\tHAVANA\tCDS\t4003\t4300\t.\t-\t1\tgene_id "ENSG002.1"; transcript_id "ENST002.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 2;
chr1\tHAVANA\tstop_codon\t4000\t4002\t.\t-\t0\tgene_id "ENSG002.1"; transcript_id "ENST002.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 2;
chr1\tHAVANA\ttranscript\t7000\t7300\t.\t+\t.\tgene_id "ENSG003.1"; transcript_id "ENST003.1"; gene_type "protein_coding"; transcript_type "protein_coding";
chr1\tHAVANA\texon\t7000\t7100\t.\t+\t.\tgene_id "ENSG003.1"; transcript_id "ENST003.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chr1\tHAVANA\texon\t7200\t7300\t.\t+\t.\tgene_id "ENSG003.1"; transcript_id "ENST003.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 2;
chr1\tHAVANA\tCDS\t7051\t7098\t.\t+\t0\tgene_id "ENSG003.1"; transcript_id "ENST003.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chr1\tHAVANA\tstop_codon\t7099\t7100\t.\t+\t0\tgene_id "ENSG003.1"; transcript_id "ENST003.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chr1\tHAVANA\tstop_codon\t7200\t7200\t.\t+\t1\tgene_id "ENSG003.1"; transcript_id "ENST003.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 2;
chr1\tHAVANA\ttranscript\t8000\t8300\t.\t+\t.\tgene_id "ENSG004.1"; transcript_id "ENST004.1"; gene_type "protein_coding"; transcript_type "protein_coding"; tag "cds_end_NF";
chr1\tHAVANA\texon\t8000\t8300\t.\t+\t.\tgene_id "ENSG004.1"; transcript_id "ENST004.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chr1\tHAVANA\tCDS\t8050\t8300\t.\t+\t0\tgene_id "ENSG004.1"; transcript_id "ENST004.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chr1\tHAVANA\ttranscript\t9000\t9300\t.\t+\t.\tgene_id "ENSG005.1"; transcript_id "ENST005.1"; gene_type "protein_coding"; transcript_type "protein_coding"; tag "cds_end_NF";
chr1\tHAVANA\texon\t9000\t9300\t.\t+\t.\tgene_id "ENSG005.1"; transcript_id "ENST005.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chr1\tHAVANA\tCDS\t9051\t9098\t.\t+\t0\tgene_id "ENSG005.1"; transcript_id "ENST005.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chr1\tHAVANA\ttranscript\t9500\t9900\t.\t-\t.\tgene_id "ENSG007.1"; transcript_id "ENST007.1"; gene_type "protein_coding"; transcript_type "protein_coding"; tag "cds_end_NF";
chr1\tHAVANA\texon\t9800\t9900\t.\t-\t.\tgene_id "ENSG007.1"; transcript_id "ENST007.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chr1\tHAVANA\texon\t9500\t9600\t.\t-\t.\tgene_id "ENSG007.1"; transcript_id "ENST007.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 2;
chr1\tHAVANA\tCDS\t9802\t9849\t.\t-\t0\tgene_id "ENSG007.1"; transcript_id "ENST007.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chrM\tENSEMBL\ttranscript\t100\t400\t.\t+\t.\tgene_id "ENSG006.1"; transcript_id "ENST006.1"; gene_type "protein_coding"; transcript_type "protein_coding";
chrM\tENSEMBL\texon\t100\t400\t.\t+\t.\tgene_id "ENSG006.1"; transcript_id "ENST006.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
chrM\tENSEMBL\tCDS\t101\t160\t.\t+\t0\tgene_id "ENSG006.1"; transcript_id "ENST006.1"; gene_type "protein_coding"; transcript_type "protein_coding"; exon_number 1;
"""

_GENCODE_GFF3 = """\
##gff-version 3
chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tID=ENSG001.1;gene_id=ENSG001.1;gene_type=protein_coding;gene_name=GENE1
chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST001.1;Parent=ENSG001.1;gene_id=ENSG001.1;transcript_id=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding
chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tID=exon:ENST001.1:1;Parent=ENST001.1;gene_id=ENSG001.1;transcript_id=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chr1\tHAVANA\texon\t1500\t2000\t.\t+\t.\tID=exon:ENST001.1:2;Parent=ENST001.1;gene_id=ENSG001.1;transcript_id=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=2
chr1\tHAVANA\tCDS\t1051\t1200\t.\t+\t0\tID=CDS:ENST001.1;Parent=ENST001.1;gene_id=ENSG001.1;transcript_id=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chr1\tHAVANA\tCDS\t1500\t2000\t.\t+\t0\tID=CDS:ENST001.1;Parent=ENST001.1;gene_id=ENSG001.1;transcript_id=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=2
chr1\tHAVANA\tstop_codon\t1998\t2000\t.\t+\t0\tID=stop_codon:ENST001.1;Parent=ENST001.1;gene_id=ENSG001.1;transcript_id=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=2
chr1\tHAVANA\tgene\t4000\t5300\t.\t-\t.\tID=ENSG002.1;gene_id=ENSG002.1;gene_type=protein_coding;gene_name=GENE2
chr1\tHAVANA\ttranscript\t4000\t5300\t.\t-\t.\tID=ENST002.1;Parent=ENSG002.1;gene_id=ENSG002.1;transcript_id=ENST002.1;gene_type=protein_coding;transcript_type=protein_coding
chr1\tHAVANA\texon\t5100\t5300\t.\t-\t.\tID=exon:ENST002.1:1;Parent=ENST002.1;gene_id=ENSG002.1;transcript_id=ENST002.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chr1\tHAVANA\texon\t4000\t4300\t.\t-\t.\tID=exon:ENST002.1:2;Parent=ENST002.1;gene_id=ENSG002.1;transcript_id=ENST002.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=2
chr1\tHAVANA\tCDS\t5100\t5299\t.\t-\t0\tID=CDS:ENST002.1;Parent=ENST002.1;gene_id=ENSG002.1;transcript_id=ENST002.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chr1\tHAVANA\tCDS\t4000\t4300\t.\t-\t1\tID=CDS:ENST002.1;Parent=ENST002.1;gene_id=ENSG002.1;transcript_id=ENST002.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=2
chr1\tHAVANA\tstop_codon\t4000\t4002\t.\t-\t0\tID=stop_codon:ENST002.1;Parent=ENST002.1;gene_id=ENSG002.1;transcript_id=ENST002.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=2
chr1\tHAVANA\ttranscript\t7000\t7300\t.\t+\t.\tID=ENST003.1;Parent=ENSG003.1;gene_id=ENSG003.1;transcript_id=ENST003.1;gene_type=protein_coding;transcript_type=protein_coding
chr1\tHAVANA\texon\t7000\t7100\t.\t+\t.\tID=exon:ENST003.1:1;Parent=ENST003.1;gene_id=ENSG003.1;transcript_id=ENST003.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chr1\tHAVANA\texon\t7200\t7300\t.\t+\t.\tID=exon:ENST003.1:2;Parent=ENST003.1;gene_id=ENSG003.1;transcript_id=ENST003.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=2
chr1\tHAVANA\tCDS\t7051\t7100\t.\t+\t0\tID=CDS:ENST003.1;Parent=ENST003.1;gene_id=ENSG003.1;transcript_id=ENST003.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chr1\tHAVANA\tCDS\t7200\t7200\t.\t+\t1\tID=CDS:ENST003.1;Parent=ENST003.1;gene_id=ENSG003.1;transcript_id=ENST003.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=2
chr1\tHAVANA\tstop_codon\t7099\t7100\t.\t+\t0\tID=stop_codon:ENST003.1:1;Parent=ENST003.1;gene_id=ENSG003.1;transcript_id=ENST003.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chr1\tHAVANA\tstop_codon\t7200\t7200\t.\t+\t1\tID=stop_codon:ENST003.1:2;Parent=ENST003.1;gene_id=ENSG003.1;transcript_id=ENST003.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=2
chr1\tHAVANA\ttranscript\t8000\t8300\t.\t+\t.\tID=ENST004.1;Parent=ENSG004.1;gene_id=ENSG004.1;transcript_id=ENST004.1;gene_type=protein_coding;transcript_type=protein_coding;tag=cds_end_NF
chr1\tHAVANA\texon\t8000\t8300\t.\t+\t.\tID=exon:ENST004.1:1;Parent=ENST004.1;gene_id=ENSG004.1;transcript_id=ENST004.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1;tag=cds_end_NF
chr1\tHAVANA\tCDS\t8050\t8300\t.\t+\t0\tID=CDS:ENST004.1;Parent=ENST004.1;gene_id=ENSG004.1;transcript_id=ENST004.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1;tag=cds_end_NF
chr1\tHAVANA\ttranscript\t9000\t9300\t.\t+\t.\tID=ENST005.1;Parent=ENSG005.1;gene_id=ENSG005.1;transcript_id=ENST005.1;gene_type=protein_coding;transcript_type=protein_coding;tag=cds_end_NF
chr1\tHAVANA\texon\t9000\t9300\t.\t+\t.\tID=exon:ENST005.1:1;Parent=ENST005.1;gene_id=ENSG005.1;transcript_id=ENST005.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1;tag=cds_end_NF
chr1\tHAVANA\tCDS\t9051\t9101\t.\t+\t0\tID=CDS:ENST005.1;Parent=ENST005.1;gene_id=ENSG005.1;transcript_id=ENST005.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1;tag=cds_end_NF
chr1\tHAVANA\ttranscript\t9500\t9900\t.\t-\t.\tID=ENST007.1;Parent=ENSG007.1;gene_id=ENSG007.1;transcript_id=ENST007.1;gene_type=protein_coding;transcript_type=protein_coding;tag=cds_end_NF
chr1\tHAVANA\texon\t9800\t9900\t.\t-\t.\tID=exon:ENST007.1:1;Parent=ENST007.1;gene_id=ENSG007.1;transcript_id=ENST007.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1;tag=cds_end_NF
chr1\tHAVANA\texon\t9500\t9600\t.\t-\t.\tID=exon:ENST007.1:2;Parent=ENST007.1;gene_id=ENSG007.1;transcript_id=ENST007.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=2;tag=cds_end_NF
chr1\tHAVANA\tCDS\t9800\t9849\t.\t-\t0\tID=CDS:ENST007.1;Parent=ENST007.1;gene_id=ENSG007.1;transcript_id=ENST007.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1;tag=cds_end_NF
chr1\tHAVANA\tCDS\t9600\t9600\t.\t-\t1\tID=CDS:ENST007.1;Parent=ENST007.1;gene_id=ENSG007.1;transcript_id=ENST007.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=2;tag=cds_end_NF
chrM\tENSEMBL\ttranscript\t100\t400\t.\t+\t.\tID=ENST006.1;Parent=ENSG006.1;gene_id=ENSG006.1;transcript_id=ENST006.1;gene_type=protein_coding;transcript_type=protein_coding
chrM\tENSEMBL\texon\t100\t400\t.\t+\t.\tID=exon:ENST006.1:1;Parent=ENST006.1;gene_id=ENSG006.1;transcript_id=ENST006.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chrM\tENSEMBL\tCDS\t101\t160\t.\t+\t0\tID=CDS:ENST006.1;Parent=ENST006.1;gene_id=ENSG006.1;transcript_id=ENST006.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
"""

_ENSEMBL_GTF = """\
chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tgene_id "ENSGE001"; gene_biotype "protein_coding";
chr1\tensembl\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "ENSGE001"; transcript_id "ENSTE001"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
chr1\tensembl\texon\t1000\t1200\t.\t+\t.\tgene_id "ENSGE001"; transcript_id "ENSTE001"; exon_number "1";
chr1\tensembl\texon\t1500\t2000\t.\t+\t.\tgene_id "ENSGE001"; transcript_id "ENSTE001"; exon_number "2";
chr1\tensembl\tCDS\t1051\t1200\t.\t+\t0\tgene_id "ENSGE001"; transcript_id "ENSTE001"; exon_number "1";
chr1\tensembl\tCDS\t1500\t1997\t.\t+\t0\tgene_id "ENSGE001"; transcript_id "ENSTE001"; exon_number "2";
chr1\tensembl\tstop_codon\t1998\t2000\t.\t+\t0\tgene_id "ENSGE001"; transcript_id "ENSTE001"; exon_number "2";
chr1\tensembl\tgene\t4000\t5300\t.\t-\t.\tgene_id "ENSGE002"; gene_biotype "protein_coding";
chr1\tensembl\ttranscript\t4000\t5300\t.\t-\t.\tgene_id "ENSGE002"; transcript_id "ENSTE002"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
chr1\tensembl\texon\t5100\t5300\t.\t-\t.\tgene_id "ENSGE002"; transcript_id "ENSTE002"; exon_number "1";
chr1\tensembl\texon\t4000\t4300\t.\t-\t.\tgene_id "ENSGE002"; transcript_id "ENSTE002"; exon_number "2";
chr1\tensembl\tCDS\t5100\t5299\t.\t-\t0\tgene_id "ENSGE002"; transcript_id "ENSTE002"; exon_number "1";
chr1\tensembl\tCDS\t4003\t4300\t.\t-\t1\tgene_id "ENSGE002"; transcript_id "ENSTE002"; exon_number "2";
chr1\tensembl\tstop_codon\t4000\t4002\t.\t-\t0\tgene_id "ENSGE002"; transcript_id "ENSTE002"; exon_number "2";
chr1\tensembl\ttranscript\t7000\t7300\t.\t+\t.\tgene_id "ENSGE003"; transcript_id "ENSTE003"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
chr1\tensembl\texon\t7000\t7100\t.\t+\t.\tgene_id "ENSGE003"; transcript_id "ENSTE003"; exon_number "1";
chr1\tensembl\texon\t7200\t7300\t.\t+\t.\tgene_id "ENSGE003"; transcript_id "ENSTE003"; exon_number "2";
chr1\tensembl\tCDS\t7051\t7098\t.\t+\t0\tgene_id "ENSGE003"; transcript_id "ENSTE003"; exon_number "1";
chr1\tensembl\tstop_codon\t7099\t7100\t.\t+\t0\tgene_id "ENSGE003"; transcript_id "ENSTE003"; exon_number "1";
chr1\tensembl\tstop_codon\t7200\t7200\t.\t+\t1\tgene_id "ENSGE003"; transcript_id "ENSTE003"; exon_number "2";
chr1\tensembl\ttranscript\t8000\t8300\t.\t+\t.\tgene_id "ENSGE004"; transcript_id "ENSTE004"; gene_biotype "protein_coding"; transcript_biotype "protein_coding"; tag "cds_end_NF";
chr1\tensembl\texon\t8000\t8300\t.\t+\t.\tgene_id "ENSGE004"; transcript_id "ENSTE004"; exon_number "1";
chr1\tensembl\tCDS\t8050\t8300\t.\t+\t0\tgene_id "ENSGE004"; transcript_id "ENSTE004"; exon_number "1";
MT\tensembl\ttranscript\t100\t400\t.\t+\t.\tgene_id "ENSGE006"; transcript_id "ENSTE006"; gene_biotype "protein_coding"; transcript_biotype "protein_coding";
MT\tensembl\texon\t100\t400\t.\t+\t.\tgene_id "ENSGE006"; transcript_id "ENSTE006"; exon_number "1";
MT\tensembl\tCDS\t101\t157\t.\t+\t0\tgene_id "ENSGE006"; transcript_id "ENSTE006"; exon_number "1";
MT\tensembl\tstop_codon\t158\t160\t.\t+\t0\tgene_id "ENSGE006"; transcript_id "ENSTE006"; exon_number "1";
"""

_ENSEMBL_GFF3 = """\
##gff-version 3
chr1\tensembl\tgene\t1000\t2000\t.\t+\t.\tID=gene:ENSGE001;biotype=protein_coding
chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:ENSTE001;Parent=gene:ENSGE001;biotype=protein_coding
chr1\tensembl\texon\t1000\t1200\t.\t+\t.\tParent=transcript:ENSTE001;rank=1;ensembl_end_phase=0
chr1\tensembl\texon\t1500\t2000\t.\t+\t.\tParent=transcript:ENSTE001;rank=2;ensembl_end_phase=0
chr1\tensembl\tCDS\t1051\t1200\t.\t+\t0\tID=CDS:ENSPE001;Parent=transcript:ENSTE001
chr1\tensembl\tCDS\t1500\t2000\t.\t+\t0\tID=CDS:ENSPE001;Parent=transcript:ENSTE001
chr1\tensembl\tgene\t4000\t5300\t.\t-\t.\tID=gene:ENSGE002;biotype=protein_coding
chr1\tensembl\tmRNA\t4000\t5300\t.\t-\t.\tID=transcript:ENSTE002;Parent=gene:ENSGE002;biotype=protein_coding
chr1\tensembl\texon\t5100\t5300\t.\t-\t.\tParent=transcript:ENSTE002;rank=1;ensembl_end_phase=2
chr1\tensembl\texon\t4000\t4300\t.\t-\t.\tParent=transcript:ENSTE002;rank=2;ensembl_end_phase=0
chr1\tensembl\tCDS\t5100\t5299\t.\t-\t0\tID=CDS:ENSPE002;Parent=transcript:ENSTE002
chr1\tensembl\tCDS\t4000\t4300\t.\t-\t1\tID=CDS:ENSPE002;Parent=transcript:ENSTE002
chr1\tensembl\tgene\t7000\t7300\t.\t+\t.\tID=gene:ENSGE003;biotype=protein_coding
chr1\tensembl\tmRNA\t7000\t7300\t.\t+\t.\tID=transcript:ENSTE003;Parent=gene:ENSGE003;biotype=protein_coding
chr1\tensembl\texon\t7000\t7100\t.\t+\t.\tParent=transcript:ENSTE003;rank=1;ensembl_end_phase=2
chr1\tensembl\texon\t7200\t7300\t.\t+\t.\tParent=transcript:ENSTE003;rank=2;ensembl_end_phase=-1
chr1\tensembl\tCDS\t7051\t7100\t.\t+\t0\tID=CDS:ENSPE003;Parent=transcript:ENSTE003
chr1\tensembl\tCDS\t7200\t7200\t.\t+\t1\tID=CDS:ENSPE003;Parent=transcript:ENSTE003
chr1\tensembl\tgene\t8000\t8300\t.\t+\t.\tID=gene:ENSGE004;biotype=protein_coding
chr1\tensembl\tmRNA\t8000\t8300\t.\t+\t.\tID=transcript:ENSTE004;Parent=gene:ENSGE004;biotype=protein_coding
chr1\tensembl\texon\t8000\t8300\t.\t+\t.\tParent=transcript:ENSTE004;rank=1;ensembl_end_phase=2
chr1\tensembl\tCDS\t8050\t8300\t.\t+\t0\tID=CDS:ENSPE004;Parent=transcript:ENSTE004
MT\tensembl\tmRNA\t100\t400\t.\t+\t.\tID=transcript:ENSTE006;Parent=gene:ENSGE006;biotype=protein_coding
MT\tensembl\texon\t100\t400\t.\t+\t.\tParent=transcript:ENSTE006;rank=1;ensembl_end_phase=-1
MT\tensembl\tCDS\t101\t160\t.\t+\t0\tID=CDS:ENSPE006;Parent=transcript:ENSTE006
"""

# Stop codons of the fixture transcripts, as 1-based start and plus strand bases: TTA is TAA on the
# minus strand (ENST002), TA + G is the split TAG (ENST003), ENST004 ends in TAA out of frame, and
# T + TA is TAA on the minus strand, split across an intron (ENST007).
_FIXTURE_BASES = {
    ("chr1", 1998): "TAA",
    ("chr1", 4000): "TTA",
    ("chr1", 7099): "TA",
    ("chr1", 7200): "G",
    ("chr1", 8298): "TAA",
    ("chr1", 9099): "TGA",
    ("chr1", 9600): "T",
    ("chr1", 9800): "TA",
    ("chrM", 158): "AGA",
    ("MT", 158): "AGA",
}


def _fasta(tmp_path, bases=None):
    """
    Writes a FASTA with chromosomes chr1, chrM, chrX, chrY and MT of 10 kb C each, with ``bases``
    ({(chromosome, 1-based start): bases}, default _FIXTURE_BASES) in it, and returns it as a
    pyfaidx.Fasta object.
    """
    seqs = {chrom: ["C"] * 10_000 for chrom in ("chr1", "chrM", "chrX", "chrY", "MT")}
    for (chrom, start), planted in (_FIXTURE_BASES if bases is None else bases).items():
        seqs[chrom][start - 1 : start - 1 + len(planted)] = list(planted)
    path = tmp_path / "genome.fa"
    path.write_text("".join(f">{chrom}\n{''.join(seq)}\n" for chrom, seq in seqs.items()))
    return Fasta(str(path))


_COMPARISON_COLUMNS = [
    "Chromosome",
    "Feature",
    "Start",
    "End",
    "Strand",
    "transcript_id",
    "gene_id",
    "exon_number",
    "has_stop_codon",
]


def _cds_exon_table(pyranges_obj):
    """Exon rows plus the coding regions (CDS rows with has_stop_codon), as extract_ptc sees them."""
    df = pyranges_obj.df
    df = df[df["Feature"].isin(["exon", "CDS"])].copy()
    df["exon_number"] = df["exon_number"].astype(int)
    df["has_stop_codon"] = df["has_stop_codon"].astype("boolean")
    # Cast away categorical dtypes so the comparison is about values, not incidental
    # category-set/order differences between how the GTF and GFF3 paths build their frames.
    for col in ["Chromosome", "Feature", "Strand", "transcript_id", "gene_id"]:
        df[col] = df[col].astype(str)
    return df[_COMPARISON_COLUMNS].sort_values(["transcript_id", "Feature", "Start"]).reset_index(drop=True)


def test_read_gff3_gencode_flavor_matches_gtf(tmp_path):
    gtf_path = _write(tmp_path, "gencode.gtf", _GENCODE_GTF)
    gff3_path = _write(tmp_path, "gencode.gff3", _GENCODE_GFF3)

    gtf_table = _cds_exon_table(nmd_scanner.scan.read_annotation(gtf_path))
    gff3_table = _cds_exon_table(nmd_scanner.scan.read_gff3(gff3_path, _fasta(tmp_path)))

    pd.testing.assert_frame_equal(gtf_table, gff3_table)


def test_read_gff3_ensembl_flavor_matches_gtf(tmp_path):
    gtf_path = _write(tmp_path, "ensembl.gtf", _ENSEMBL_GTF)
    gff3_path = _write(tmp_path, "ensembl.gff3", _ENSEMBL_GFF3)

    gtf_table = _cds_exon_table(nmd_scanner.scan.read_annotation(gtf_path))
    gff3_table = _cds_exon_table(nmd_scanner.scan.read_gff3(gff3_path, _fasta(tmp_path)))

    pd.testing.assert_frame_equal(gtf_table, gff3_table)


def test_read_gff3_stop_codon_from_sequence(tmp_path):
    """
    An Ensembl GFF3 transcript has a stop codon if the last 3 bases of its CDS are a stop codon, in
    frame or not, unless its last coding exon ends mid-codon (ensembl_end_phase 1 or 2). On MT, the
    vertebrate mitochondrial code applies: AGA is a stop codon, TGA is not. The CDS stays as it is.
    """
    # chromosome, transcript start, transcript, CDS length, ensembl_end_phase, last 3 CDS bases
    cases = [
        ("chr1", 100, "TAA", 60, -1, "TAA"),
        ("chr1", 300, "AGA", 60, -1, "AGA"),
        ("chr1", 500, "MID_CODON", 60, 2, "TAA"),
        ("chr1", 700, "OUT_OF_FRAME", 61, -1, "TGA"),
        ("MT", 100, "MT_AGA", 60, -1, "AGA"),
        ("MT", 300, "MT_TGA", 60, -1, "TGA"),
    ]
    rows, bases = [], {}
    for chrom, start, tx, length, end_phase, codon in cases:
        rows += [
            f"{chrom}\tensembl\tmRNA\t{start}\t{start + 100}\t.\t+\t.\tID=transcript:{tx};Parent=gene:G{tx};biotype=protein_coding",
            f"{chrom}\tensembl\texon\t{start}\t{start + 100}\t.\t+\t.\tParent=transcript:{tx};rank=1;ensembl_end_phase={end_phase}",
            f"{chrom}\tensembl\tCDS\t{start + 1}\t{start + length}\t.\t+\t0\tID=CDS:P{tx};Parent=transcript:{tx}",
        ]
        bases[(chrom, start + length - 2)] = codon
    gff3_path = _write(tmp_path, "ensembl.gff3", "##gff-version 3\n" + "\n".join(rows) + "\n")

    df = nmd_scanner.scan.read_gff3(gff3_path, _fasta(tmp_path, bases)).df

    assert set(df["Feature"]) == {"exon", "CDS"}
    cds = df[df["Feature"] == "CDS"]
    assert sorted(zip(cds["transcript_id"], cds["Start"], cds["End"])) == sorted(
        (tx, start, start + length) for _, start, tx, length, _, _ in cases
    )
    assert df.loc[df["Feature"] == "exon", "has_stop_codon"].isna().all()
    has_stop = dict(zip(cds["transcript_id"], cds["has_stop_codon"]))
    assert has_stop == {
        "TAA": True,
        "AGA": False,
        "MID_CODON": False,
        "OUT_OF_FRAME": True,
        "MT_AGA": True,
        "MT_TGA": False,
    }


def test_read_gff3_gives_coding_regions_without_stop_codon_rows(tmp_path):
    fasta = _fasta(tmp_path)
    for name, content, expected in [
        (
            "gencode.gff3",
            _GENCODE_GFF3,
            {"ENST001.1", "ENST002.1", "ENST003.1"},
        ),
        ("ensembl.gff3", _ENSEMBL_GFF3, {"ENSTE001", "ENSTE002", "ENSTE003", "ENSTE006"}),
    ]:
        df = nmd_scanner.scan.read_gff3(_write(tmp_path, name, content), fasta).df
        assert set(df["Feature"]) == {"exon", "CDS"}
        assert df.loc[df["Feature"] == "exon", "has_stop_codon"].isna().all()
        cds = df[df["Feature"] == "CDS"]
        assert set(cds.loc[cds["has_stop_codon"].astype(bool), "transcript_id"]) == expected


def test_read_annotation_reassigns_gff3_exon_numbers(tmp_path):
    without_numbers = re.sub(r";exon_number=\d+", "", _GENCODE_GFF3)
    assert "exon_number" not in without_numbers
    fasta = _fasta(tmp_path)
    reassigned = nmd_scanner.scan.read_annotation(
        _write(tmp_path, "a.gff3", without_numbers), fasta, reassign_exons=True
    )
    gtf = nmd_scanner.scan.read_annotation(_write(tmp_path, "a.gtf", _GENCODE_GTF))
    pd.testing.assert_frame_equal(_cds_exon_table(reassigned), _cds_exon_table(gtf))


@pytest.mark.parametrize(
    ("gtf", "gff3", "chrom_m", "n_results"),
    [(_GENCODE_GTF, _GENCODE_GFF3, "chrM", 9), (_ENSEMBL_GTF, _ENSEMBL_GFF3, "MT", 7)],
    ids=["gencode", "ensembl"],
)
@pytest.mark.parametrize("reassign_exons", [False, True])
def test_main_gives_the_same_results_for_gtf_and_gff3(tmp_path, gtf, gff3, chrom_m, n_results, reassign_exons):
    """Variants in every fixture transcript, e.g. the split stop codon and the cds_end_NF ones."""
    _fasta(tmp_path)
    variants = [
        ("chr1", 1100, "CC", "C"),
        ("chr1", 1600, "C", "A"),
        ("chr1", 4100, "C", "T"),
        ("chr1", 5200, "C", "A"),
        ("chr1", 7060, "C", "A"),
        ("chr1", 8100, "C", "T"),
        ("chr1", 9060, "C", "A"),
        ("chr1", 9820, "C", "A"),
        (chrom_m, 120, "C", "A"),
    ]
    vcf = _write(
        tmp_path,
        "variants.vcf",
        "".join(
            f"{chrom}\t{pos}\tv{i}\t{ref}\t{alt}\t.\tPASS\t.\n" for i, (chrom, pos, ref, alt) in enumerate(variants)
        ),
    )
    fasta = str(tmp_path / "genome.fa")
    out = str(tmp_path / "out.csv")
    via_gtf = main(vcf, _write(tmp_path, "a.gtf", gtf), fasta, out, reassign_exons=reassign_exons)
    via_gff3 = main(
        vcf, None, fasta, out, reassign_exons=reassign_exons, annotation_path=_write(tmp_path, "a.gff3", gff3)
    )
    # one row per variant: the Ensembl fixture has no transcript at 9060 and 9820
    assert len(via_gtf) == n_results
    pd.testing.assert_frame_equal(via_gtf, via_gff3)


def test_read_annotation_gff3_needs_fasta(tmp_path):
    gff3_path = _write(tmp_path, "ensembl.gff3", _ENSEMBL_GFF3)
    with pytest.raises(ValueError, match="FASTA"):
        nmd_scanner.scan.read_annotation(gff3_path)


def test_read_gff3_drops_id_and_parent_columns(tmp_path):
    """
    GFF3's own ID/Parent columns must not leak into the returned frame: downstream code joins
    the CDS table against the VCF table (which has its own unrelated "ID" column, the variant
    ID), and a leftover GFF3 "ID" column would silently take priority in that join and replace
    the variant ID with a GFF3 feature ID such as "CDS:ENST00000643195.1".
    """
    gencode_path = _write(tmp_path, "gencode.gff3", _GENCODE_GFF3)
    ensembl_path = _write(tmp_path, "ensembl.gff3", _ENSEMBL_GFF3)

    fasta = _fasta(tmp_path)
    for path in (gencode_path, ensembl_path):
        df = nmd_scanner.scan.read_gff3(path, fasta).df
        assert "ID" not in df.columns
        assert "Parent" not in df.columns


def test_read_annotation_reads_gzipped_gff3(tmp_path):
    gz_path = tmp_path / "gencode.gff3.gz"
    with gzip.open(gz_path, "wt") as fh:
        fh.write(_GENCODE_GFF3)
    plain_path = _write(tmp_path, "gencode.gff3", _GENCODE_GFF3)

    fasta = _fasta(tmp_path)
    gz_table = _cds_exon_table(nmd_scanner.scan.read_annotation(str(gz_path), fasta))
    plain_table = _cds_exon_table(nmd_scanner.scan.read_gff3(plain_path, fasta))
    pd.testing.assert_frame_equal(gz_table, plain_table)


def test_read_gff3_gencode_resolves_par_y_like_ids_via_hierarchy(tmp_path):
    """
    GENCODE GFF3 chrY PAR transcripts reuse the chrX transcript_id/gene_id attribute value;
    only ID/Parent (which carry a "_PAR_Y"-like suffix) disambiguate them. Two transcripts here
    share the flat transcript_id/gene_id "ENST001.1"/"ENSG001.1", distinguished only by ID/Parent.
    """
    content = """\
##gff-version 3
chrY\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tID=ENSG001.1_PAR_Y;gene_id=ENSG001.1;gene_type=protein_coding;gene_name=GENE1
chrY\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST001.1_PAR_Y;Parent=ENSG001.1_PAR_Y;gene_id=ENSG001.1;transcript_id=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding
chrY\tHAVANA\texon\t1000\t1200\t.\t+\t.\tID=exon:ENST001.1_PAR_Y:1;Parent=ENST001.1_PAR_Y;gene_id=ENSG001.1;transcript_id=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chrY\tHAVANA\tCDS\t1050\t1200\t.\t+\t0\tID=CDS:ENST001.1_PAR_Y;Parent=ENST001.1_PAR_Y;gene_id=ENSG001.1;transcript_id=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chrX\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tID=ENSG001.1;gene_id=ENSG001.1;gene_type=protein_coding;gene_name=GENE1
chrX\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST001.1;Parent=ENSG001.1;gene_id=ENSG001.1;transcript_id=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding
chrX\tHAVANA\texon\t1000\t1200\t.\t+\t.\tID=exon:ENST001.1:1;Parent=ENST001.1;gene_id=ENSG001.1;transcript_id=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chrX\tHAVANA\tCDS\t1050\t1200\t.\t+\t0\tID=CDS:ENST001.1;Parent=ENST001.1;gene_id=ENSG001.1;transcript_id=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
"""
    gff3_path = _write(tmp_path, "par_y.gff3", content)

    df = nmd_scanner.scan.read_gff3(gff3_path, _fasta(tmp_path)).df
    cds = df[df["Feature"] == "CDS"]

    transcript_ids = set(cds["transcript_id"])
    assert transcript_ids == {"ENST001.1_PAR_Y", "ENST001.1"}

    gene_ids = dict(zip(cds["transcript_id"], cds["gene_id"]))
    assert gene_ids == {"ENST001.1_PAR_Y": "ENSG001.1_PAR_Y", "ENST001.1": "ENSG001.1"}

    # each copy keeps its own CDS, not merged or cross-contaminated with the other
    assert len(cds) == 2
    assert dict(zip(cds["transcript_id"], cds["Chromosome"].astype(str))) == {
        "ENST001.1_PAR_Y": "chrY",
        "ENST001.1": "chrX",
    }


def test_read_gff3_unrecognized_flavor_raises(tmp_path):
    content = """\
##gff-version 3
chr1\tsource\tgene\t1000\t2000\t.\t+\t.\tID=G1
chr1\tsource\ttranscript\t1000\t2000\t.\t+\t.\tID=T1;Parent=G1
"""
    gff3_path = _write(tmp_path, "unknown.gff3", content)
    with pytest.raises(ValueError, match="Unrecognized GFF3 flavor"):
        nmd_scanner.scan.read_gff3(gff3_path, _fasta(tmp_path))


# GFF3 exon numbers, ids, file names


def test_compute_exon_numbers_equal_the_gtf_numbers_on_both_strands(gtf_path):
    """The numbers computed from genomic order are the ones the chr18 GTF carries."""
    gtf = pr.read_gtf(gtf_path)
    # without the GTF's own numbers, so that the test sees only what compute_exon_numbers computes
    computed = nmd_scanner.compute_exon_numbers(pr.PyRanges(gtf.df.drop(columns="exon_number"))).df
    expected = pd.to_numeric(gtf.df["exon_number"], errors="coerce").astype("Int64")
    selected = computed["Feature"].isin(["exon", "CDS", "stop_codon"])
    assert set(computed.loc[selected, "Strand"]) == {"+", "-"}
    pd.testing.assert_series_equal(computed.loc[selected, "exon_number"], expected[selected])


def test_read_gff3_ensembl_takes_exon_numbers_from_rank_not_from_compute_exon_numbers(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise AssertionError("compute_exon_numbers must not run when the exons have a rank")

    monkeypatch.setattr(nmd_scanner.scan, "compute_exon_numbers", fail)
    gff3_path = _write(tmp_path, "ensembl.gff3", _ENSEMBL_GFF3)
    df = nmd_scanner.scan.read_gff3(gff3_path, _fasta(tmp_path)).df
    numbers = df.groupby(["transcript_id", "Feature"])["exon_number"].apply(list)
    assert numbers[("ENSTE002", "exon")] == [1, 2]
    assert numbers[("ENSTE002", "CDS")] == [1, 2]
    assert numbers[("ENSTE003", "CDS")] == [1, 2]


def test_read_gff3_ensembl_exon_number_is_the_rank_attribute(tmp_path):
    """The rank is used as given, also where it differs from the genomic order."""
    content = """\
##gff-version 3
chr1\tensembl\tmRNA\t1000\t2000\t.\t+\t.\tID=transcript:T1;Parent=gene:G1;biotype=protein_coding
chr1\tensembl\texon\t1000\t1200\t.\t+\t.\tParent=transcript:T1;rank=7
chr1\tensembl\texon\t1500\t2000\t.\t+\t.\tParent=transcript:T1;rank=8
chr1\tensembl\tCDS\t1051\t1200\t.\t+\t0\tID=CDS:P1;Parent=transcript:T1
chr1\tensembl\tCDS\t1500\t1900\t.\t+\t0\tID=CDS:P1;Parent=transcript:T1
"""
    df = nmd_scanner.scan.read_gff3(_write(tmp_path, "rank.gff3", content), _fasta(tmp_path)).df
    assert df[df["Feature"] == "CDS"]["exon_number"].tolist() == [7, 8]
    assert df[df["Feature"] == "exon"]["exon_number"].tolist() == [7, 8]


def test_read_gff3_ensembl_without_rank_computes_exon_numbers(tmp_path):
    content = _ENSEMBL_GFF3.replace("rank=1;", "").replace("rank=2;", "")
    assert "rank=" not in content
    fasta = _fasta(tmp_path)
    without_rank = nmd_scanner.scan.read_gff3(_write(tmp_path, "no_rank.gff3", content), fasta)
    with_rank = nmd_scanner.scan.read_gff3(_write(tmp_path, "rank.gff3", _ENSEMBL_GFF3), fasta)
    pd.testing.assert_frame_equal(_cds_exon_table(without_rank), _cds_exon_table(with_rank))


def test_compute_exon_numbers_cds_takes_the_exon_with_the_most_overlap():
    df = pd.DataFrame(
        {
            "Chromosome": ["chr1"] * 4,
            "Start": [100, 200, 190, 195],
            "End": [200, 300, 205, 240],
            "Strand": ["-"] * 4,
            "Feature": ["exon", "exon", "CDS", "CDS"],
            "transcript_id": ["T1"] * 4,
        }
    )
    out = nmd_scanner.compute_exon_numbers(pr.PyRanges(df)).df.sort_values("Start")
    assert out.loc[out["Feature"] == "exon", "exon_number"].tolist() == [2, 1]
    # CDS 190-205 overlaps exon 100-200 by 10 and exon 200-300 by 5. CDS 195-240 overlaps exon
    # 100-200 by 5 and exon 200-300 by 40, so the first overlapping exon is not the answer there.
    assert out.loc[out["Feature"] == "CDS", "exon_number"].tolist() == [2, 1]


def test_detect_annotation_format_accepts_path_objects():
    assert nmd_scanner.scan.detect_annotation_format(Path("a.GTF")) == "gtf"
    assert nmd_scanner.scan.detect_annotation_format(Path("a.gff3.gz")) == "gff3"


def test_read_annotation_accepts_path_objects(tmp_path):
    gtf_path = _write(tmp_path, "a.gtf", _GENCODE_GTF)
    gff3_path = _write(tmp_path, "a.gff3", _GENCODE_GFF3)
    fasta = _fasta(tmp_path)
    gtf_table = _cds_exon_table(nmd_scanner.scan.read_annotation(Path(gtf_path), fasta))
    gff3_table = _cds_exon_table(nmd_scanner.scan.read_annotation(Path(gff3_path), fasta))
    pd.testing.assert_frame_equal(gtf_table, gff3_table)


def test_read_annotation_fmt_overrides_the_file_suffix(tmp_path):
    plain_path = _write(tmp_path, "plain.txt", _GENCODE_GTF)
    with pytest.raises(ValueError, match="Cannot detect annotation format"):
        nmd_scanner.scan.read_annotation(plain_path)
    via_fmt = nmd_scanner.scan.read_annotation(plain_path, fmt="gtf").df
    pd.testing.assert_frame_equal(via_fmt, nmd_scanner.scan.read_annotation(_write(tmp_path, "a.gtf", _GENCODE_GTF)).df)
    with pytest.raises(ValueError, match="Unknown annotation format"):
        nmd_scanner.scan.read_annotation(plain_path, fmt="bed")


def test_read_gff3_gencode_keeps_the_transcript_id_and_gene_id_attributes(tmp_path):
    """GENCODE lift37 ids carry a _N suffix on the attributes, but not in ID/Parent; the GTF has the suffix."""
    content = """\
##gff-version 3
chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tID=ENSG001.1;gene_id=ENSG001.1_9;gene_type=protein_coding
chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST001.1;Parent=ENSG001.1;gene_id=ENSG001.1_9;transcript_id=ENST001.1_2;gene_type=protein_coding;transcript_type=protein_coding
chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tID=exon:ENST001.1:1;Parent=ENST001.1;gene_id=ENSG001.1_9;transcript_id=ENST001.1_2;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chr1\tHAVANA\tCDS\t1050\t1200\t.\t+\t0\tID=CDS:ENST001.1;Parent=ENST001.1;gene_id=ENSG001.1_9;transcript_id=ENST001.1_2;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
"""
    df = nmd_scanner.scan.read_gff3(_write(tmp_path, "lift37.gff3", content), _fasta(tmp_path)).df
    assert set(df["transcript_id"]) == {"ENST001.1_2"}
    assert set(df["gene_id"]) == {"ENSG001.1_9"}


def test_read_gff3_gencode_fills_missing_id_attributes_from_the_hierarchy(tmp_path):
    content = """\
##gff-version 3
chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tID=ENSG001.1;gene_type=protein_coding
chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST001.1;Parent=ENSG001.1;gene_type=protein_coding;transcript_type=protein_coding
chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tID=exon:ENST001.1:1;Parent=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
"""
    df = nmd_scanner.scan.read_gff3(_write(tmp_path, "no_ids.gff3", content), _fasta(tmp_path)).df
    assert df["transcript_id"].tolist() == ["ENST001.1"]
    assert df["gene_id"].tolist() == ["ENSG001.1"]


def test_read_gff3_gencode_fills_id_attributes_missing_on_some_rows(tmp_path):
    """ENST000.1 has the attributes, ENST001.1 has none: the hierarchy fills in only ENST001.1."""
    content = """\
##gff-version 3
chr1\tHAVANA\tgene\t100\t500\t.\t+\t.\tID=ENSG000.1;gene_id=ENSG000.1_5;gene_type=protein_coding
chr1\tHAVANA\ttranscript\t100\t500\t.\t+\t.\tID=ENST000.1;Parent=ENSG000.1;gene_id=ENSG000.1_5;transcript_id=ENST000.1_5;gene_type=protein_coding;transcript_type=protein_coding
chr1\tHAVANA\texon\t100\t500\t.\t+\t.\tID=exon:ENST000.1:1;Parent=ENST000.1;gene_id=ENSG000.1_5;transcript_id=ENST000.1_5;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tID=ENSG001.1;gene_type=protein_coding
chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST001.1;Parent=ENSG001.1;gene_type=protein_coding;transcript_type=protein_coding
chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tID=exon:ENST001.1:1;Parent=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
"""
    df = nmd_scanner.scan.read_gff3(_write(tmp_path, "some_ids.gff3", content), _fasta(tmp_path)).df
    assert df["transcript_id"].tolist() == ["ENST000.1_5", "ENST001.1"]
    assert df["gene_id"].tolist() == ["ENSG000.1_5", "ENSG001.1"]


def test_a_gtf_with_a_gff3_name_raises_an_error_naming_the_file(tmp_path):
    path = _write(tmp_path, "really_a_gtf.gff3", _GENCODE_GTF)
    with pytest.raises(ValueError, match=r"really_a_gtf\.gff3.*as GFF3.*GTF"):
        nmd_scanner.scan.read_annotation(path, _fasta(tmp_path))


def test_a_gff3_with_a_gtf_name_raises_an_error_naming_the_file(tmp_path):
    path = _write(tmp_path, "really_a_gff3.gtf", _GENCODE_GFF3)
    with pytest.raises(ValueError, match=r"really_a_gff3\.gtf.*as GTF.*GFF3"):
        nmd_scanner.scan.read_annotation(path)
