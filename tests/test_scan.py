# Import dependencies
import gzip
import logging
import random
import re
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from pyfaidx import Fasta

import nmd_scanner
from nmd_scanner.cli import main

# pytest-fixtures as inputs for the tests

VCF_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


@pytest.fixture(scope="session")
def gff3_path():
    return "resources/chr18.gff3.gz"


@pytest.fixture(scope="session")
def vcf_path():
    return "resources/part-00241-61a0abbf-fbf9-444f-8287-4e46ad4b9b7b-c000.vcf"


@pytest.fixture(scope="session")
def fasta_path():
    return "resources/chr18.fa.gz"


# Create the test functions


# Test reading VCF file
def test_read_vcf_file(vcf_path):

    df = nmd_scanner.scan.read_vcf(vcf_path)
    assert isinstance(df, pd.DataFrame)
    assert df.shape[0] > 0
    assert list(df.columns) == ["Chromosome", "Start", "End", "ID", "Ref", "Alt"]
    pd.testing.assert_index_equal(df.index, pd.RangeIndex(len(df)))


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


@pytest.mark.parametrize("alt", ["G,<*>", "G]chr2:1],T", "<INS:ME|ALU>,T"])
def test_read_vcf_rejects_multiallelic_records_with_symbolic_alleles_and_breakends(tmp_path, alt):
    vcf = tmp_path / "multiallelic.vcf"
    vcf.write_text(VCF_HEADER + f"chr1\t100\tv1\tA\t{alt}\t.\t.\t.\n")
    with pytest.raises(ValueError, match="1 multi-allelic record"):
        nmd_scanner.scan.read_vcf(str(vcf))


# VCF 4.3 allows "|" in the ID of a symbolic allele and in a contig name, which a breakend names
@pytest.mark.parametrize("alt", ["<INS:ME|ALU>", "G]gi|123|:100]", "[gi|123|:100[G"])
def test_read_vcf_accepts_a_symbolic_allele_or_breakend_with_a_pipe(tmp_path, alt):
    vcf = tmp_path / "pipe.vcf"
    vcf.write_text(VCF_HEADER + f"chr1\t100\tv1\tG\t{alt}\t.\t.\t.\n")
    df = nmd_scanner.scan.read_vcf(str(vcf))
    assert df["Alt"].tolist() == [alt]


def test_read_vcf_accepts_single_allelic(tmp_path):
    vcf = tmp_path / "single.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\tv1\tA\tT\t.\t.\t.\n"
        "chr1\t200\tv2\tC\tG\t.\t.\t.\n"
    )
    df = nmd_scanner.scan.read_vcf(str(vcf))
    assert df.shape[0] == 2


# A column with any non-numeric value stays text even without dtype=str, so the all-numeric ID pair
# ("007", "0123") is what pins the int inference; ("12345", "NA") pins the NA parsing.
@pytest.mark.parametrize("ids", [("007", "0123"), ("12345", "NA")])
def test_read_vcf_keeps_text_fields_as_written(tmp_path, ids):
    vcf = tmp_path / "text_fields.vcf"
    vcf.write_text(VCF_HEADER + f"01\t100\t{ids[0]}\tA\tNA\t.\tPASS\tNA\n01\t200\t{ids[1]}\tC\tG\t50\t.\t.\n")
    df = nmd_scanner.scan.read_vcf(str(vcf))
    assert (df["Chromosome"] == "01").all()
    assert df["ID"].tolist() == list(ids)
    assert df["Alt"].tolist() == ["NA", "G"]
    assert df["Start"].tolist() == [99, 199]
    assert df["End"].tolist() == [100, 200]


def test_read_vcf_keeps_a_dot_in_id_and_alt(tmp_path):
    """polars-bio gives "" for "."; read_vcf turns it back into "."."""
    vcf = tmp_path / "dots.vcf"
    vcf.write_text(VCF_HEADER + "chr1\t100\t.\tA\t.\t.\t.\t.\nchr1\t200\tv2\tC\tG\t.\t.\t.\n")
    df = nmd_scanner.scan.read_vcf(str(vcf))
    assert df["ID"].tolist() == [".", "v2"]
    assert df["Alt"].tolist() == [".", "G"]


def test_read_vcf_drops_qual_filter_and_info(tmp_path):
    vcf = tmp_path / "all_fields.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\tv1\tA\tT\t50\tPASS\tDP=12\n"
    )
    df = nmd_scanner.scan.read_vcf(str(vcf))
    assert list(df.columns) == ["Chromosome", "Start", "End", "ID", "Ref", "Alt"]


def test_read_vcf_end_comes_from_ref_not_from_info_end(tmp_path):
    """polars-bio takes the end of a symbolic allele from INFO END; read_vcf keeps Start plus the REF length."""
    vcf = tmp_path / "symbolic.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t100\tdel\tN\t<DEL>\t.\t.\tEND=200\n"
        "chr1\t300\tindel\tACG\tA\t.\t.\t.\n"
    )
    df = nmd_scanner.scan.read_vcf(str(vcf))
    assert df["Start"].tolist() == [99, 299]
    assert df["End"].tolist() == [100, 302]


@pytest.mark.parametrize(
    "header",
    ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n", "##fileformat=VCFv4.2\n", ""],
    ids=["without_fileformat_line", "without_chrom_line", "without_header"],
)
def test_read_vcf_without_header_raises_an_error_naming_the_file(tmp_path, header):
    vcf = tmp_path / "no_header.vcf"
    vcf.write_text(header + "chr1\t100\tv1\tA\tT\t.\t.\t.\n")
    with pytest.raises(ValueError, match=r"no_header\.vcf.*needs its header.*##fileformat.*#CHROM"):
        nmd_scanner.scan.read_vcf(str(vcf))


def test_read_vcf_gives_the_header_hint_only_for_a_missing_header(tmp_path):
    vcf = tmp_path / "trailing_blank.vcf"
    vcf.write_text(VCF_HEADER + "chr1\t100\tv1\tA\tT\t.\t.\t.\n\n")
    with pytest.raises(ValueError, match=r"trailing_blank\.vcf") as error:
        nmd_scanner.scan.read_vcf(str(vcf))
    assert "needs its header" not in str(error.value)


def test_read_vcf_of_a_directory_raises_is_a_directory_error(tmp_path):
    directory = tmp_path / "variants.vcf"
    directory.mkdir()
    with pytest.raises(IsADirectoryError):
        nmd_scanner.scan.read_vcf(str(directory))


def test_read_vcf_reads_plain_and_gzip_files_alike(tmp_path):
    content = VCF_HEADER + "chr1\t100\tv1\tA\tT\t.\t.\t.\nchr1\t200\tv2\tCA\tC\t.\t.\t.\n"
    plain = tmp_path / "variants.vcf"
    plain.write_text(content)
    gz = tmp_path / "variants.vcf.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(content)
    expected = nmd_scanner.scan.read_vcf(str(plain))
    assert len(expected) == 2
    pd.testing.assert_frame_equal(nmd_scanner.scan.read_vcf(str(gz)), expected)


def test_read_vcf_missing_file_raises_file_not_found(tmp_path):
    with pytest.raises(FileNotFoundError, match="missing.vcf"):
        nmd_scanner.scan.read_vcf(str(tmp_path / "missing.vcf"))


def _coding_sequence(coding, fasta):
    coding = coding.sort_values("Start")
    seq = "".join(fasta[c][s:e].seq.upper() for c, s, e in zip(coding["Chromosome"], coding["Start"], coding["End"]))
    return str(Seq(seq).reverse_complement()) if coding["Strand"].iloc[0] == "-" else seq


def test_read_annotation_gives_whole_coding_regions_on_real_transcripts(gff3_path, fasta_path):
    fasta = Fasta(fasta_path)
    annotation = nmd_scanner.scan.read_annotation(gff3_path, fasta)
    cds = annotation[annotation["Feature"] == "CDS"]

    # a stop codon split across an intron: its last base (+ strand) or its last 2 bases (- strand) make up the
    # CDS row of the last coding exon
    for transcript_id, strand, last_row in [("ENST00000399496.8", "+", -1), ("ENST00000454642.3", "-", 0)]:
        coding = cds[cds["transcript_id"] == transcript_id].sort_values("Start")
        assert set(coding["Strand"]) == {strand}
        assert coding["End"].iloc[last_row] - coding["Start"].iloc[last_row] < 3
        assert coding["has_stop_codon"].all()
        seq = _coding_sequence(coding, fasta)
        assert len(seq) % 3 == 0
        assert seq[-3:] in {"TAA", "TAG", "TGA"}

    # a cds_end_NF transcript has no stop codon
    cds_end_nf = cds["tag"].fillna("").str.contains("cds_end_NF")
    assert cds_end_nf.any()
    assert not cds.loc[cds_end_nf, "has_stop_codon"].any()


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
    out1 = nmd_scanner.compute_exon_numbers(df1)

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

    out2 = nmd_scanner.compute_exon_numbers(df2)

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
    out4 = nmd_scanner.compute_exon_numbers(df4)

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
    out3 = nmd_scanner.compute_exon_numbers(df3)

    ex_txA = out3[(out3.Feature == "exon") & (out3.transcript_id == "TXA")].sort_values("Start")
    ex_txB = out3[(out3.Feature == "exon") & (out3.transcript_id == "TXB")].sort_values("Start")
    assert list(ex_txA["exon_number"]) == [1, 2]
    assert list(ex_txB["exon_number"]) == [1, 2, 3]

    # TODO: maybe add the edge case if:
    # 1. CDS does not overlap any exon --> should not crash but exon_number should stay missing
    # 2. CDS overlaps two exons --> should it inherit the exon_number with the maximum overlap??


def test_compute_exon_numbers_takes_a_dataframe_with_any_index():
    df = pd.DataFrame(
        [
            ["chr1", 100, 200, "+", "exon", "TX1"],
            ["chr1", 300, 400, "+", "exon", "TX1"],
            ["chr1", 320, 400, "+", "CDS", "TX1"],
        ],
        columns=["Chromosome", "Start", "End", "Strand", "Feature", "transcript_id"],
        index=[7, 7, 7],
    )
    out = nmd_scanner.compute_exon_numbers(df)
    assert isinstance(out, pd.DataFrame)
    pd.testing.assert_index_equal(out.index, pd.RangeIndex(3))
    assert out["exon_number"].tolist() == [1, 2, 2]


def test_compute_exon_numbers_with_str_exon_number_column():
    """
    A caller can give exon_number as ``str``, with missing values on features that have no exon
    number. Computed exon numbers must be ints, not written into the str column.
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
    out = nmd_scanner.compute_exon_numbers(df)

    exons = out[out.Feature == "exon"].sort_values("Start")
    assert list(exons["exon_number"]) == [1, 2]
    cds = out[out.Feature == "CDS"].iloc[0]
    assert cds["exon_number"] == 2
    assert not isinstance(cds["exon_number"], str)
    assert pd.isna(out[out.Feature == "gene"].iloc[0]["exon_number"])
    assert pd.api.types.is_integer_dtype(out["exon_number"])


# Test annotation format detection and dispatch


def test_detect_annotation_format():
    assert nmd_scanner.scan.detect_annotation_format("annotation.gff3") == "gff3"
    assert nmd_scanner.scan.detect_annotation_format("annotation.gff3.gz") == "gff3"
    assert nmd_scanner.scan.detect_annotation_format("annotation.gff") == "gff3"
    assert nmd_scanner.scan.detect_annotation_format("annotation.gff.gz") == "gff3"
    assert nmd_scanner.scan.detect_annotation_format("ANNOTATION.GFF3") == "gff3"

    with pytest.raises(ValueError, match="Cannot detect annotation format"):
        nmd_scanner.scan.detect_annotation_format("annotation.txt")
    for name in ("annotation.gtf", "annotation.gtf.gz", "ANNOTATION.GTF"):
        with pytest.raises(ValueError, match=f"Cannot read '{name}': GTF input is no longer supported"):
            nmd_scanner.scan.detect_annotation_format(name)


def _write(tmp_path, name, content):
    path = tmp_path / name
    path.write_text(content)
    return str(path)


# A tiny GENCODE-style GFF3, and the equivalent Ensembl-style GFF3, for the same loci: a + and a -
# strand transcript, a stop codon split across an intron (ENST003.1), a transcript without stop
# codon (cds_end_NF, ENST004.1) whose last 3 bases read TAA out of frame, and a chrM transcript
# whose CDS ends in AGA (ENST006.1). The GFF3 CDS includes the stop codon. The bases are in
# _FIXTURE_BASES. _GENCODE_ROWS and _ENSEMBL_ROWS hold the rows that read_gff3 gives for them:
# the coding regions and has_stop_codon of the GTF of the same loci.
# Ensembl: the GFF3 has no stop_codon rows, read_gff3 takes has_stop_codon from the FASTA. AGA is a
# stop codon on MT. ENSTE004 has ensembl_end_phase 2, so it has no stop codon.
# GENCODE: there is no stop codon for AGA on chrM, since the GFF3 has no stop_codon rows there.
# ENST005.1 (+ strand) and ENST007.1 (- strand, split across an intron) are tagged cds_end_NF, but
# their CDS ends in a complete stop codon without stop_codon rows; the GENCODE GTF has these 3 bases
# as UTR, and read_gff3 removes them from the CDS.

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


def _cds_exon_table(annotation):
    """Exon rows plus the coding regions (CDS rows with has_stop_codon), as extract_ptc sees them."""
    df = annotation[annotation["Feature"].isin(["exon", "CDS"])].copy()
    df["exon_number"] = df["exon_number"].astype(int)
    df["has_stop_codon"] = df["has_stop_codon"].astype("boolean")
    # Cast away categorical dtypes so the comparison is about values, not incidental
    # category-set/order differences.
    for col in ["Chromosome", "Feature", "Strand", "transcript_id", "gene_id"]:
        df[col] = df[col].astype(str)
    return df[_COMPARISON_COLUMNS].sort_values(["transcript_id", "Feature", "Start"]).reset_index(drop=True)


def _rows(annotation):
    """
    The rows of _cds_exon_table as (transcript_id, Feature, Start, End, exon_number, has_stop_codon)
    tuples, has_stop_codon None on the exon rows, and {transcript_id: (gene_id, Chromosome, Strand)}.
    """
    table = _cds_exon_table(annotation)
    rows = [
        (tx, feature, start, end, number, None if pd.isna(has_stop) else bool(has_stop))
        for tx, feature, start, end, number, has_stop in table[
            ["transcript_id", "Feature", "Start", "End", "exon_number", "has_stop_codon"]
        ].itertuples(index=False)
    ]
    transcripts = dict(zip(table["transcript_id"], zip(table["gene_id"], table["Chromosome"], table["Strand"])))
    return rows, transcripts


# 0-based half-open; a coding region is the CDS plus the stop codon bases of its exon
_GENCODE_ROWS = [
    ("ENST001.1", "CDS", 1050, 1200, 1, True),
    ("ENST001.1", "CDS", 1499, 2000, 2, True),
    ("ENST001.1", "exon", 999, 1200, 1, None),
    ("ENST001.1", "exon", 1499, 2000, 2, None),
    ("ENST002.1", "CDS", 3999, 4300, 2, True),
    ("ENST002.1", "CDS", 5099, 5299, 1, True),
    ("ENST002.1", "exon", 3999, 4300, 2, None),
    ("ENST002.1", "exon", 5099, 5300, 1, None),
    ("ENST003.1", "CDS", 7050, 7100, 1, True),
    # second base of the split stop codon, in an exon without other coding bases
    ("ENST003.1", "CDS", 7199, 7200, 2, True),
    ("ENST003.1", "exon", 6999, 7100, 1, None),
    ("ENST003.1", "exon", 7199, 7300, 2, None),
    ("ENST004.1", "CDS", 8049, 8300, 1, False),
    ("ENST004.1", "exon", 7999, 8300, 1, None),
    # cds_end_NF: the CDS without its last 3 bases, which read TGA
    ("ENST005.1", "CDS", 9050, 9098, 1, False),
    ("ENST005.1", "exon", 8999, 9300, 1, None),
    ("ENST006.1", "CDS", 100, 160, 1, False),
    ("ENST006.1", "exon", 99, 400, 1, None),
    # cds_end_NF: without the stop codon split across the intron, so without the CDS row of exon 2
    ("ENST007.1", "CDS", 9801, 9849, 1, False),
    ("ENST007.1", "exon", 9499, 9600, 2, None),
    ("ENST007.1", "exon", 9799, 9900, 1, None),
]
_GENCODE_TRANSCRIPTS = {
    "ENST001.1": ("ENSG001.1", "chr1", "+"),
    "ENST002.1": ("ENSG002.1", "chr1", "-"),
    "ENST003.1": ("ENSG003.1", "chr1", "+"),
    "ENST004.1": ("ENSG004.1", "chr1", "+"),
    "ENST005.1": ("ENSG005.1", "chr1", "+"),
    "ENST006.1": ("ENSG006.1", "chrM", "+"),
    "ENST007.1": ("ENSG007.1", "chr1", "-"),
}
_ENSEMBL_ROWS = [
    ("ENSTE001", "CDS", 1050, 1200, 1, True),
    ("ENSTE001", "CDS", 1499, 2000, 2, True),
    ("ENSTE001", "exon", 999, 1200, 1, None),
    ("ENSTE001", "exon", 1499, 2000, 2, None),
    ("ENSTE002", "CDS", 3999, 4300, 2, True),
    ("ENSTE002", "CDS", 5099, 5299, 1, True),
    ("ENSTE002", "exon", 3999, 4300, 2, None),
    ("ENSTE002", "exon", 5099, 5300, 1, None),
    ("ENSTE003", "CDS", 7050, 7100, 1, True),
    ("ENSTE003", "CDS", 7199, 7200, 2, True),
    ("ENSTE003", "exon", 6999, 7100, 1, None),
    ("ENSTE003", "exon", 7199, 7300, 2, None),
    ("ENSTE004", "CDS", 8049, 8300, 1, False),
    ("ENSTE004", "exon", 7999, 8300, 1, None),
    ("ENSTE006", "CDS", 100, 160, 1, True),
    ("ENSTE006", "exon", 99, 400, 1, None),
]
_ENSEMBL_TRANSCRIPTS = {
    "ENSTE001": ("ENSGE001", "chr1", "+"),
    "ENSTE002": ("ENSGE002", "chr1", "-"),
    "ENSTE003": ("ENSGE003", "chr1", "+"),
    "ENSTE004": ("ENSGE004", "chr1", "+"),
    "ENSTE006": ("ENSGE006", "MT", "+"),
}


@pytest.mark.parametrize(
    ("gff3", "expected_rows", "expected_transcripts"),
    [(_GENCODE_GFF3, _GENCODE_ROWS, _GENCODE_TRANSCRIPTS), (_ENSEMBL_GFF3, _ENSEMBL_ROWS, _ENSEMBL_TRANSCRIPTS)],
    ids=["gencode", "ensembl"],
)
def test_read_gff3_gives_the_coding_regions_and_exons_of_the_fixture(
    tmp_path, gff3, expected_rows, expected_transcripts
):
    annotation = nmd_scanner.scan.read_gff3(_write(tmp_path, "a.gff3", gff3), _fasta(tmp_path))
    assert isinstance(annotation, pd.DataFrame)
    pd.testing.assert_index_equal(annotation.index, pd.RangeIndex(len(annotation)))
    assert annotation["exon_number"].dtype == "Int64"
    assert _rows(annotation) == (expected_rows, expected_transcripts)


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

    df = nmd_scanner.scan.read_gff3(gff3_path, _fasta(tmp_path, bases))

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
        df = nmd_scanner.scan.read_gff3(_write(tmp_path, name, content), fasta)
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
    assert _rows(reassigned) == (_GENCODE_ROWS, _GENCODE_TRANSCRIPTS)


# variant_id, transcript_id, has_stop_codon, ref_cds_len, alt_cds_len, ref_cds_info: the values that the GTF of
# the same loci gave
_GENCODE_RESULTS = [
    ("v0", "ENST001.1", True, 651, 650, [(1, 150), (2, 501)]),
    ("v1", "ENST001.1", True, 651, 651, [(1, 150), (2, 501)]),
    ("v2", "ENST002.1", True, 501, 501, [(1, 200), (2, 301)]),
    ("v3", "ENST002.1", True, 501, 501, [(1, 200), (2, 301)]),
    ("v4", "ENST003.1", True, 51, 51, [(1, 50), (2, 1)]),
    ("v5", "ENST004.1", False, 251, 251, [(1, 251)]),
    ("v6", "ENST005.1", False, 48, 48, [(1, 48)]),
    ("v8", "ENST006.1", False, 60, 60, [(1, 60)]),
    ("v7", "ENST007.1", False, 48, 48, [(1, 48)]),
]
_ENSEMBL_RESULTS = [
    ("v0", "ENSTE001", True, 651, 650, [(1, 150), (2, 501)]),
    ("v1", "ENSTE001", True, 651, 651, [(1, 150), (2, 501)]),
    ("v2", "ENSTE002", True, 501, 501, [(1, 200), (2, 301)]),
    ("v3", "ENSTE002", True, 501, 501, [(1, 200), (2, 301)]),
    ("v4", "ENSTE003", True, 51, 51, [(1, 50), (2, 1)]),
    ("v5", "ENSTE004", False, 251, 251, [(1, 251)]),
    # the Ensembl fixture has no transcript at 9060 and 9820, and AGA is a stop codon on MT
    ("v8", "ENSTE006", True, 60, 60, [(1, 60)]),
]


@pytest.mark.parametrize(
    ("gff3", "chrom_m", "expected"),
    [(_GENCODE_GFF3, "chrM", _GENCODE_RESULTS), (_ENSEMBL_GFF3, "MT", _ENSEMBL_RESULTS)],
    ids=["gencode", "ensembl"],
)
def test_main_gives_the_results_of_the_fixture(tmp_path, gff3, chrom_m, expected):
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
        VCF_HEADER
        + "".join(
            f"{chrom}\t{pos}\tv{i}\t{ref}\t{alt}\t.\tPASS\t.\n" for i, (chrom, pos, ref, alt) in enumerate(variants)
        ),
    )
    fasta = str(tmp_path / "genome.fa")
    gff3_path = _write(tmp_path, "a.gff3", gff3)
    results = main(vcf, gff3_path, fasta, str(tmp_path / "out.csv"))
    columns = ["variant_id", "transcript_id", "has_stop_codon", "ref_cds_len", "alt_cds_len", "ref_cds_info"]
    assert [tuple(row) for row in results[columns].itertuples(index=False)] == expected
    reassigned = main(vcf, gff3_path, fasta, str(tmp_path / "reassigned.csv"), reassign_exons=True)
    pd.testing.assert_frame_equal(reassigned, results)


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
        df = nmd_scanner.scan.read_gff3(path, fasta)
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

    df = nmd_scanner.scan.read_gff3(gff3_path, _fasta(tmp_path))
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


def test_read_gff3_percent_decodes_attribute_values(tmp_path):
    """GFF3 escapes ";", "=", "&" and "," in attribute values as %3B, %3D, %26 and %2C."""
    content = """\
##gff-version 3
chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST001.1;Parent=ENSG001.1;gene_id=ENSG%3D001%26x;transcript_id=ENST%3B001%2C1;gene_type=protein_coding;transcript_type=protein_coding
chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tID=exon:ENST001.1:1;Parent=ENST001.1;gene_id=ENSG%3D001%26x;transcript_id=ENST%3B001%2C1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chr1\tHAVANA\tCDS\t1051\t1200\t.\t+\t0\tID=CDS:ENST001.1;Parent=ENST001.1;gene_id=ENSG%3D001%26x;transcript_id=ENST%3B001%2C1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
"""
    df = nmd_scanner.scan.read_gff3(_write(tmp_path, "escaped.gff3", content), _fasta(tmp_path))
    assert set(df["transcript_id"]) == {"ENST;001,1"}
    assert set(df["gene_id"]) == {"ENSG=001&x"}


def test_read_gff3_unrecognized_flavor_raises(tmp_path):
    content = """\
##gff-version 3
chr1\tsource\tgene\t1000\t2000\t.\t+\t.\tID=G1
chr1\tsource\ttranscript\t1000\t2000\t.\t+\t.\tID=T1;Parent=G1
"""
    gff3_path = _write(tmp_path, "unknown.gff3", content)
    with pytest.raises(ValueError, match="Unrecognized GFF3 flavor"):
        nmd_scanner.scan.read_gff3(gff3_path, _fasta(tmp_path))


_CDS_LINE = "chr1\tHAVANA\tCDS\t1051\t1200\t.\t+\t0\tID=CDS:ENST001.1;"


# polars-bio skips these lines without an error
@pytest.mark.parametrize(
    "line",
    [
        _CDS_LINE.replace("\t.\t+", "\tabc\t+"),
        _CDS_LINE.replace("1051", "abc"),
        _CDS_LINE.replace("1051", "-1051"),
        _CDS_LINE.replace("1051", str(2**32)),
        _CDS_LINE.replace("\t", " "),
    ],
    ids=["text_score", "text_start", "negative_start", "start_of_2_to_the_32", "space_separated"],
)
def test_read_gff3_raises_for_a_line_that_polars_bio_skips(tmp_path, line):
    assert _GENCODE_GFF3.count(_CDS_LINE) == 1
    data_lines = len(_GENCODE_GFF3.splitlines()) - 1
    gff3_path = _write(tmp_path, "malformed.gff3", _GENCODE_GFF3.replace(_CDS_LINE, line))
    with pytest.raises(
        ValueError, match=rf"malformed\.gff3.*polars-bio read {data_lines - 1} of its {data_lines} data lines"
    ):
        nmd_scanner.scan.read_gff3(gff3_path, _fasta(tmp_path))


def test_read_gff3_does_not_count_blank_lines_comments_and_the_fasta_section(tmp_path):
    lines = _GENCODE_GFF3.splitlines()
    content = (
        "\n".join(lines[:3] + ["", "  \t", "# a comment", "###"] + lines[3:] + ["##FASTA", ">chr1", "ACGT"]) + "\n\n"
    )
    fasta = _fasta(tmp_path)
    expected = nmd_scanner.scan.read_gff3(_write(tmp_path, "plain.gff3", _GENCODE_GFF3), fasta)
    pd.testing.assert_frame_equal(nmd_scanner.scan.read_gff3(_write(tmp_path, "extra.gff3", content), fasta), expected)


@pytest.mark.parametrize(
    ("start", "end", "shown"),
    [("1201", "1200", "chr1:1201-1200"), ("0", "1200", "chr1:4294967296-1200")],
    ids=["start_after_end", "start_of_0"],
)
def test_read_gff3_raises_for_a_start_after_the_end(tmp_path, start, end, shown):
    line = _CDS_LINE.replace("\t1051\t1200\t", f"\t{start}\t{end}\t")
    gff3_path = _write(tmp_path, "reversed.gff3", _GENCODE_GFF3.replace(_CDS_LINE, line))
    with pytest.raises(
        ValueError, match=rf"reversed\.gff3.*1 row\(s\) have a start after their end, e.g. the CDS row at {shown}"
    ):
        nmd_scanner.scan.read_gff3(gff3_path, _fasta(tmp_path))


def test_read_gff3_warns_if_a_gencode_gff3_has_no_stop_codon_rows(tmp_path, caplog):
    without_stop_codons = "".join(line for line in _GENCODE_GFF3.splitlines(True) if "\tstop_codon\t" not in line)
    fasta = _fasta(tmp_path)
    with caplog.at_level(logging.WARNING, logger="nmd_scanner.scan"):
        annotation = nmd_scanner.scan.read_gff3(_write(tmp_path, "no_stop.gff3", without_stop_codons), fasta)
    assert "No stop_codon rows found next to the CDS rows" in caplog.text
    assert not annotation["has_stop_codon"].any()

    caplog.clear()
    with caplog.at_level(logging.WARNING, logger="nmd_scanner.scan"):
        nmd_scanner.scan.read_gff3(_write(tmp_path, "stop.gff3", _GENCODE_GFF3), fasta)
    assert "No stop_codon rows" not in caplog.text


def test_read_gff3_raises_an_error_naming_the_file_if_polars_bio_cannot_read_it(tmp_path):
    path = tmp_path / "random.gff3"
    path.write_bytes(random.Random(0).randbytes(300))
    with pytest.raises(ValueError, match=r"Cannot read '.*random\.gff3' as GFF3 \("):
        nmd_scanner.scan.read_gff3(str(path), _fasta(tmp_path))


# GFF3 exon numbers, ids, file names


def test_compute_exon_numbers_equal_the_annotated_numbers_on_both_strands(gff3_path, fasta_path):
    """The numbers computed from genomic order are the ones the chr18 GFF3 carries."""
    annotation = nmd_scanner.scan.read_gff3(gff3_path, Fasta(fasta_path))
    # without the annotated numbers, so that the test sees only what compute_exon_numbers computes
    computed = nmd_scanner.compute_exon_numbers(annotation.drop(columns="exon_number"))
    expected = annotation["exon_number"].astype("Int64")
    assert set(computed["Feature"]) == {"exon", "CDS"}
    assert set(computed["Strand"]) == {"+", "-"}
    pd.testing.assert_series_equal(computed["exon_number"], expected)


def test_read_gff3_ensembl_takes_exon_numbers_from_rank_not_from_compute_exon_numbers(tmp_path, monkeypatch):
    def fail(*args, **kwargs):
        raise AssertionError("compute_exon_numbers must not run when the exons have a rank")

    monkeypatch.setattr(nmd_scanner.scan, "compute_exon_numbers", fail)
    gff3_path = _write(tmp_path, "ensembl.gff3", _ENSEMBL_GFF3)
    df = nmd_scanner.scan.read_gff3(gff3_path, _fasta(tmp_path))
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
    df = nmd_scanner.scan.read_gff3(_write(tmp_path, "rank.gff3", content), _fasta(tmp_path))
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
    out = nmd_scanner.compute_exon_numbers(df).sort_values("Start")
    assert out.loc[out["Feature"] == "exon", "exon_number"].tolist() == [2, 1]
    # CDS 190-205 overlaps exon 100-200 by 10 and exon 200-300 by 5. CDS 195-240 overlaps exon
    # 100-200 by 5 and exon 200-300 by 40, so the first overlapping exon is not the answer there.
    assert out.loc[out["Feature"] == "CDS", "exon_number"].tolist() == [2, 1]


def test_detect_annotation_format_accepts_path_objects():
    assert nmd_scanner.scan.detect_annotation_format(Path("a.gff3.gz")) == "gff3"
    with pytest.raises(ValueError, match="GTF input is no longer supported"):
        nmd_scanner.scan.detect_annotation_format(Path("a.GTF"))


def test_read_annotation_accepts_path_objects(tmp_path):
    gff3_path = _write(tmp_path, "a.gff3", _GENCODE_GFF3)
    fasta = _fasta(tmp_path)
    assert _rows(nmd_scanner.scan.read_annotation(Path(gff3_path), fasta)) == (_GENCODE_ROWS, _GENCODE_TRANSCRIPTS)


def test_read_annotation_fmt_overrides_the_file_suffix(tmp_path):
    plain_path = _write(tmp_path, "plain.txt", _GENCODE_GFF3)
    fasta = _fasta(tmp_path)
    with pytest.raises(ValueError, match="Cannot detect annotation format"):
        nmd_scanner.scan.read_annotation(plain_path, fasta)
    via_fmt = nmd_scanner.scan.read_annotation(plain_path, fasta, fmt="gff3")
    pd.testing.assert_frame_equal(
        via_fmt, nmd_scanner.scan.read_annotation(_write(tmp_path, "a.gff3", _GENCODE_GFF3), fasta)
    )
    with pytest.raises(ValueError, match="Unknown annotation format 'bed', expected 'gff3'"):
        nmd_scanner.scan.read_annotation(plain_path, fasta, fmt="bed")


def test_read_gff3_gencode_keeps_the_transcript_id_and_gene_id_attributes(tmp_path):
    """GENCODE lift37 ids carry a _N suffix on the attributes, but not in ID/Parent; the GTF has the suffix."""
    content = """\
##gff-version 3
chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tID=ENSG001.1;gene_id=ENSG001.1_9;gene_type=protein_coding
chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST001.1;Parent=ENSG001.1;gene_id=ENSG001.1_9;transcript_id=ENST001.1_2;gene_type=protein_coding;transcript_type=protein_coding
chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tID=exon:ENST001.1:1;Parent=ENST001.1;gene_id=ENSG001.1_9;transcript_id=ENST001.1_2;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
chr1\tHAVANA\tCDS\t1050\t1200\t.\t+\t0\tID=CDS:ENST001.1;Parent=ENST001.1;gene_id=ENSG001.1_9;transcript_id=ENST001.1_2;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
"""
    df = nmd_scanner.scan.read_gff3(_write(tmp_path, "lift37.gff3", content), _fasta(tmp_path))
    assert set(df["transcript_id"]) == {"ENST001.1_2"}
    assert set(df["gene_id"]) == {"ENSG001.1_9"}


def test_read_gff3_gencode_fills_missing_id_attributes_from_the_hierarchy(tmp_path):
    content = """\
##gff-version 3
chr1\tHAVANA\tgene\t1000\t2000\t.\t+\t.\tID=ENSG001.1;gene_type=protein_coding
chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST001.1;Parent=ENSG001.1;gene_type=protein_coding;transcript_type=protein_coding
chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tID=exon:ENST001.1:1;Parent=ENST001.1;gene_type=protein_coding;transcript_type=protein_coding;exon_number=1
"""
    df = nmd_scanner.scan.read_gff3(_write(tmp_path, "no_ids.gff3", content), _fasta(tmp_path))
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
    df = nmd_scanner.scan.read_gff3(_write(tmp_path, "some_ids.gff3", content), _fasta(tmp_path))
    assert df["transcript_id"].tolist() == ["ENST000.1_5", "ENST001.1"]
    assert df["gene_id"].tolist() == ["ENSG000.1_5", "ENSG001.1"]


_GTF_ROWS = """\
chr1\tHAVANA\ttranscript\t1000\t2000\t.\t+\t.\tgene_id "ENSG001.1"; transcript_id "ENST001.1"; gene_type "protein_coding";
chr1\tHAVANA\texon\t1000\t1200\t.\t+\t.\tgene_id "ENSG001.1"; transcript_id "ENST001.1"; exon_number 1;
chr1\tHAVANA\tCDS\t1051\t1200\t.\t+\t0\tgene_id "ENSG001.1"; transcript_id "ENST001.1"; exon_number 1;
"""


def test_a_gtf_with_a_gff3_name_raises_an_error_naming_the_file(tmp_path):
    path = _write(tmp_path, "really_a_gtf.gff3", _GTF_ROWS)
    with pytest.raises(
        ValueError, match=r"really_a_gtf\.gff3.*as GFF3.*If it is a GTF: GTF input is no longer supported"
    ):
        nmd_scanner.scan.read_annotation(path, _fasta(tmp_path))


@pytest.mark.parametrize("name", ["a.gtf", "a.gtf.gz", "A.GTF"])
def test_read_annotation_rejects_a_gtf_file_name(tmp_path, name):
    """The file name decides, although this file holds a GFF3."""
    path = _write(tmp_path, name, _GENCODE_GFF3)
    with pytest.raises(ValueError) as error:
        nmd_scanner.scan.read_annotation(path, _fasta(tmp_path))
    assert str(error.value) == (
        f"Cannot read {path!r}: GTF input is no longer supported. Use the GFF3 of the same GENCODE or Ensembl release."
    )


def test_read_annotation_rejects_fmt_gtf(tmp_path):
    path = _write(tmp_path, "a.gff3", _GENCODE_GFF3)
    with pytest.raises(ValueError, match=r"with fmt='gtf': GTF input is no longer supported"):
        nmd_scanner.scan.read_annotation(path, _fasta(tmp_path), fmt="gtf")
