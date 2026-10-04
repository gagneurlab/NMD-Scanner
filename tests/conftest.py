"""Fixtures shared by the test modules."""

import pytest


@pytest.fixture
def intergenic_vcf(tmp_path):
    """A VCF with one variant on chr18 outside every CDS."""

    path = tmp_path / "intergenic.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr18\t1000\tintergenic\tA\tT\t.\t.\t.\n"
    )
    return str(path)


@pytest.fixture
def reference_mismatch_vcf(tmp_path):
    """A VCF with one variant inside a CDS whose REF does not match the FASTA (it is ATG there)."""

    path = tmp_path / "mismatch.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr18\t21383518\tmismatch\tCCC\tTTT\t.\t.\t.\n"
    )
    return str(path)
