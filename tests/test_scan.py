# Import dependencies
import os
import tempfile
import pytest
import nmd_scanner
import gzip

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

    keys = list(fasta.keys())  # convert to list for safety
    assert isinstance(keys, list)
    assert len(keys) > 0