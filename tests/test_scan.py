# Import dependencies
import pandas as pd
import pyranges as pr
import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

import nmd_scanner
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
