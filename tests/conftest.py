"""
Pytest configuration and fixtures for SDRF Refiner tests.
"""

import pytest
from pathlib import Path
import pandas as pd
from dataclasses import dataclass
from typing import Optional, Dict, Any


# Path to test fixtures
FIXTURES_DIR = Path(__file__).parent / "fixtures"


@pytest.fixture
def fixtures_dir():
    """Return path to test fixtures directory."""
    return FIXTURES_DIR


@pytest.fixture
def sample_sdrf_path(fixtures_dir):
    """Return path to sample SDRF file."""
    return fixtures_dir / "sample_sdrf.tsv"


@pytest.fixture
def minimal_sdrf_path(fixtures_dir):
    """Return path to minimal SDRF file."""
    return fixtures_dir / "minimal_sdrf.tsv"


@pytest.fixture
def sample_sdrf_df():
    """Return a sample SDRF DataFrame for testing."""
    return pd.DataFrame({
        "source name": ["sample1", "sample2", "sample3"],
        "comment[data file]": ["sample1.raw", "sample2.raw", "sample3.raw"],
        "comment[instrument]": [
            "NT=Q Exactive HF;AC=MS:1002523",
            "NT=Q Exactive HF;AC=MS:1002523",
            "NT=Q Exactive HF;AC=MS:1002523",
        ],
        "comment[dissociation method]": [
            "NT=HCD;AC=MS:1000422",
            "NT=HCD;AC=MS:1000422",
            "NT=HCD;AC=MS:1000422",
        ],
        "comment[precursor mass tolerance]": ["10 ppm", "10 ppm", "10 ppm"],
        "comment[fragment mass tolerance]": ["0.02 Da", "0.02 Da", "0.02 Da"],
    })


@dataclass
class MockAnalysisResult:
    """Mock AnalysisResult for testing confidence scoring."""
    instrument_model: Dict[str, Any]
    fragmentation_type: str
    precursor_tolerance_ppm: Optional[float]
    fragment_tolerance_ppm: Optional[float] = None
    fragment_tolerance_da: Optional[float] = None


@pytest.fixture
def mock_analysis_result():
    """Return a mock analysis result."""
    return MockAnalysisResult(
        instrument_model={"name": "Q Exactive HF", "accession": "MS:1002523"},
        fragmentation_type="HR_HCD",
        precursor_tolerance_ppm=10.0,
        fragment_tolerance_ppm=15.0,
    )


@pytest.fixture
def mock_analysis_results():
    """Return a list of mock analysis results."""
    return [
        MockAnalysisResult(
            instrument_model={"name": "Q Exactive HF", "accession": "MS:1002523"},
            fragmentation_type="HR_HCD",
            precursor_tolerance_ppm=10.0,
            fragment_tolerance_ppm=15.0,
        ),
        MockAnalysisResult(
            instrument_model={"name": "Q Exactive HF", "accession": "MS:1002523"},
            fragmentation_type="HR_HCD",
            precursor_tolerance_ppm=10.5,
            fragment_tolerance_ppm=14.8,
        ),
        MockAnalysisResult(
            instrument_model={"name": "Q Exactive HF", "accession": "MS:1002523"},
            fragmentation_type="HR_HCD",
            precursor_tolerance_ppm=9.8,
            fragment_tolerance_ppm=15.2,
        ),
    ]


# Add pytest mark for slow tests
def pytest_addoption(parser):
    """Add custom pytest options."""
    parser.addoption(
        "--runslow",
        action="store_true",
        default=False,
        help="run slow tests",
    )


def pytest_configure(config):
    """Configure pytest markers."""
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    """Skip slow tests unless --runslow is specified."""
    if config.getoption("--runslow"):
        return

    skip_slow = pytest.mark.skip(reason="need --runslow option to run")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
