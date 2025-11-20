# Import dependencies
import pytest
import nmd_scanner
import pandas as pd
import pyranges as pr

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


# Test reading GTF file
def test_read_gtf_file(gtf_path):
    gr = nmd_scanner.scan.read_gtf(gtf_path)
    assert gr is not None
    assert gr.df.shape[0] > 0
    assert "Chromosome" in gr.df.columns
    print(gr.df.head())


# Test reading FASTA file

def test_read_fasta_file(fasta_path):
    fasta = nmd_scanner.scan.read_fasta(fasta_path)
    assert fasta is not None

    keys = list(fasta.keys())
    assert isinstance(keys, list)
    assert len(keys) > 0

def compute_exon_numbers():

    # On + Strand: Smallest exon number is the Start, Largest exon number is the End.
    df1 = pd.DataFrame([
        ["chr1", 100, 200, "+", "exon", "TX1", "G1"],  # exon → 1
        ["chr1", 300, 400, "+", "exon", "TX1", "G1"],  # exon → 2
        ["chr1", 320, 400, "+", "CDS", "TX1", "G1"],  # CDS on exon 2
    ], columns=["Chromosome", "Start", "End", "Strand", "Feature", "transcript_id", "gene_id"])
    out1 = compute_exon_numbers(pr.PyRanges(df1)).df

    tx1_exons = out1[(out1.Feature == "exon") & (out1.transcript_id == "TX1")].sort_values("Start")
    assert list(tx1_exons["exon_number"]) == [1, 2]
    tx1_cds = out1[(out1.Feature == "CDS") & (out1.transcript_id == "TX1")].iloc[0]
    assert tx1_cds["exon_number"] == 2


    # On - Strand: Smallest exon number is the Start, Largest exon number is the End.
    df2 = pd.DataFrame([
        ["chr1", 100, 200, "-", "exon", "TX2", "G2"],  # exon_number → 2 (reverse order)
        ["chr1", 300, 400, "-", "exon", "TX2", "G2"],  # exon_number → 1
        ["chr1", 120, 180, "-", "CDS", "TX2", "G2"], # CDS on exon 2
    ], columns=["Chromosome", "Start", "End", "Strand", "Feature", "transcript_id", "gene_id"])

    out2 = compute_exon_numbers(pr.PyRanges(df2)).df

    tx2_exons = out2[(out2.Feature == "exon") & (out2.transcript_id == "TX2")].sort_values("Start")
    assert list(tx2_exons["exon_number"]) == [2, 1]
    tx2_cds = out2[(out2.Feature == "CDS") & (out2.transcript_id == "TX2")].iloc[0]
    assert tx2_cds["exon_number"] == 2


    # multiple CDS sequences on minus strand
    df4 = pd.DataFrame([
        # exons
        ["chr1", 30, 50, "-", "exon", "TX2b", "G2b"],  # exon number 4
        ["chr1", 100, 200, "-", "exon", "TX2b", "G2b"],  # exon number 3
        ["chr1", 300, 350, "-", "exon", "TX2b", "G2b"],  # exon number 2
        ["chr1", 500, 600, "-", "exon", "TX2b", "G2b"],  # exon number 1
        # cds segments
        ["chr1", 500, 590, "-", "CDS", "TX2b", "G2b"],  # exon_number 1
        ["chr1", 300, 350, "-", "CDS", "TX2b", "G2b"],  # exon_number 2
        ["chr1", 130, 200, "-", "CDS", "TX2b", "G2b"],  # exon_number 3
    ], columns=["Chromosome", "Start", "End", "Strand", "Feature", "transcript_id", "gene_id"])
    out4 = compute_exon_numbers(pr.PyRanges(df4)).df

    # check exons exon-numbers
    tx2b_exons = out4[(out4.Feature == "exon") & (out4.transcript_id == "TX2b")].sort_values("Start")
    assert list(tx2b_exons["exon_number"]) == [4, 3, 2, 1]

    # check CDS exon-numbers
    assert int(out4[(out4.Feature == "CDS") & (out4.transcript_id == "TX2b") &
                     (out4.Start == 520) & (out4.End == 590)]["exon_number"].iloc[0]) == 1
    assert int(out4[(out4.Feature == "CDS") & (out4.transcript_id == "TX2b") &
                     (out4.Start == 310) & (out4.End == 340)]["exon_number"].iloc[0]) == 2
    assert int(out4[(out4.Feature == "CDS") & (out4.transcript_id == "TX2b") &
                     (out4.Start == 130) & (out4.End == 180)]["exon_number"].iloc[0]) == 3
    assert int(out4[(out4.Feature == "CDS") & (out4.transcript_id == "TX2b") &
                     (out4.Start == 270) & (out4.End == 330)]["exon_number"].iloc[0]) == 2


    # two different transripts (should be numbered independently)
    df3 = pd.DataFrame([
        # TXA (+)
        ["chr1", 100, 150, "+", "exon", "TXA", "GA"],  # → 1
        ["chr1", 200, 250, "+", "exon", "TXA", "GA"],  # → 2
        # TXB (+)
        ["chr1", 500, 600, "+", "exon", "TXB", "GB"],  # → 1
        ["chr1", 700, 800, "+", "exon", "TXB", "GB"],  # → 2
        ["chr1", 900, 1000, "+", "exon", "TXB", "GB"],  # → 3
    ], columns=["Chromosome", "Start", "End", "Strand", "Feature", "transcript_id", "gene_id"])
    out3 = compute_exon_numbers(pr.PyRanges(df3)).df

    ex_txA = out3[(out3.Feature == "exon") & (out3.transcript_id == "TXA")].sort_values("Start")
    ex_txB = out3[(out3.Feature == "exon") & (out3.transcript_id == "TXB")].sort_values("Start")
    assert list(ex_txA["exon_number"]) == [1, 2]
    assert list(ex_txB["exon_number"]) == [1, 2, 3]

    # TODO: maybe add the edge case if:
        # 1. CDS does not overlap any exon --> should not crash but exon_number should stay missing
        # 2. CDS overlaps two exons --> should it inherit the exon_number with the maximum overlap??


