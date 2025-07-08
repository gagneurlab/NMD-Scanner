# Import dependencies
import os
import tempfile
import pytest
import nmd_scanner

# Create the test functions

# Tests for reading in the VCF file
def test_list_vcf_files():
    vcf_dir = "data/distinct_variants.valid_snp_indel.vcf"  # Directory to the VCF files
    files = nmd_scanner.scan.list_vcf_files(vcf_dir)
    assert len(files) > 0

def test_list_vcf_files_empty():
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            nmd_scanner.scan.list_vcf_files(tmpdir)
            assert False, "Expected FileNotFoundError"
        except FileNotFoundError:
            pass  # Expected

def test_single_vcf_file():
    with tempfile.NamedTemporaryFile(suffix=".vcf") as tmp:
        files, is_single = nmd_scanner.scan.list_vcf_files(tmp.name)
        assert is_single
        assert len(files) == 1
        assert files[0] == tmp.name

def test_directory_with_multiple_vcf_files():
    with tempfile.TemporaryDirectory() as tmpdir:
        paths = []
        for i in range(3):
            path = nmd_scanner.scan.os.path.join(tmpdir, f"file_{i}.vcf")
            with open(path, "w") as f:
                f.write("##dummy header\n")
            paths.append(path)

        files, is_single = nmd_scanner.scan.list_vcf_files(tmpdir)
        assert not is_single
        assert len(files) == 3
        assert set(files) == set(paths)

def test_directory_with_no_vcf_files():
    with tempfile.TemporaryDirectory() as tmpdir:
        path = nmd_scanner.scan.os.path.join(tmpdir, "not_a_vcf.txt")
        with open(path, "w") as f:
            f.write("Not a VCF")
        with pytest.raises(FileNotFoundError):
            nmd_scanner.scan.list_vcf_files(tmpdir)

def test_file_is_not_vcf():
    with tempfile.NamedTemporaryFile(suffix=".txt") as tmp:
        with pytest.raises(ValueError):
            nmd_scanner.scan.list_vcf_files(tmp.name)

def test_path_does_not_exist():
    with pytest.raises(FileNotFoundError):
        nmd_scanner.scan.list_vcf_files("/non/existing/path/to/file.vcf")

def test_directory_with_mixed_files():
    with tempfile.TemporaryDirectory() as tmpdir:
        vcf_path = nmd_scanner.scan.os.path.join(tmpdir, "valid.vcf")
        txt_path = nmd_scanner.scan.os.path.join(tmpdir, "note.txt")
        with open(vcf_path, "w") as f:
            f.write("##VCF content\n")
        with open(txt_path, "w") as f:
            f.write("Not a VCF")

        files, is_single = nmd_scanner.scan.list_vcf_files(tmpdir)
        assert not is_single
        assert len(files) == 1
        assert files[0] == vcf_path


# Test reading VCF file

def test_read_real_vcf_file():
    # Adjust this path to point to your actual file location
    real_vcf_path = "data/distinct_variants.valid_snp_indel.vcf/part-00241-61a0abbf-fbf9-444f-8287-4e46ad4b9b7b-c000.vcf"

    gr = nmd_scanner.scan.read_vcf(real_vcf_path)

    # Basic assertions to make sure the file is parsed correctly
    assert gr is not None
    assert gr.df.shape[0] > 0
    assert "Chromosome" in gr.df.columns
    assert "Start" in gr.df.columns
    assert "End" in gr.df.columns

    print(gr.df.head())
    print(gr.df.shape)


# Test reading GTF file

def test_read_gtf_file():
    gtf_path = "/s/genomes/Gencode/Gencode_human/release_40/gencode.v40.annotation.gtf.gz"

    gr = nmd_scanner.scan.read_gtf(gtf_path)

    assert gr is not None
    assert "Chromosome" in gr.df.columns
    assert gr.df.shape[0] > 0

    print(gr.df.head())