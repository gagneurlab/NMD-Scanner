import pytest

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