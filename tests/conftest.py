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
    """
    A VCF with one variant whose REF CCC does not match the FASTA, which has CAT at chr18:21383518-21383520.
    The variant overlaps the first 2 bases of the GREB1L CDS, which starts at 21383519.
    """

    path = tmp_path / "mismatch.vcf"
    path.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr18\t21383518\tmismatch\tCCC\tTTT\t.\t.\t.\n"
    )
    return str(path)
