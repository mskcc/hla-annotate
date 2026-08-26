from pathlib import Path
import pytest

FIXTURES = Path(__file__).parent / "fixtures"

@pytest.fixture
def fixtures_dir():
    return FIXTURES

@pytest.fixture
def pgroup_file(fixtures_dir):
    return fixtures_dir / "hla_nom_p_mini.txt"

@pytest.fixture
def sample1_est(fixtures_dir):
    return fixtures_dir / "sample1_A.est.txt"

@pytest.fixture
def sample2_est(fixtures_dir):
    return fixtures_dir / "sample2_A.est.txt"

@pytest.fixture
def sample3_est(fixtures_dir):
    return fixtures_dir / "sample3_A.est.txt"

@pytest.fixture
def sample4_est(fixtures_dir):
    return fixtures_dir / "sample4_A.est.txt"

@pytest.fixture
def sample5_est(fixtures_dir):
    return fixtures_dir / "sample5_A.est.txt"

@pytest.fixture
def sample3_final(fixtures_dir):
    return fixtures_dir / "sample3_final.result.txt"

@pytest.fixture
def samplename2_final(fixtures_dir):
    return fixtures_dir / "samplename2_final.result.txt"
