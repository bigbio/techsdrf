"""
Tests for confidence scoring utilities.
"""

import pytest
from unittest.mock import MagicMock
from dataclasses import dataclass
from typing import Optional, Dict, Any

from sdrf_refiner.report.confidence import (
    calculate_confidence,
    get_confidence_factors,
    _instrument_confidence,
    _fragmentation_confidence,
    _tolerance_confidence,
    _score_to_label,
)


@dataclass
class MockAnalysisResult:
    """Mock AnalysisResult for testing."""
    instrument_model: Dict[str, Any]
    fragmentation_type: str
    precursor_tolerance_ppm: Optional[float]
    fragment_tolerance_ppm: Optional[float] = None
    fragment_tolerance_da: Optional[float] = None


class TestCalculateConfidence:
    """Tests for calculate_confidence function."""

    def test_empty_results(self):
        """Test confidence with empty results."""
        score = calculate_confidence([], "instrument")
        assert score == 0.0

    def test_single_file_base_confidence(self):
        """Test base confidence with single file."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
            ),
        ]

        # Single file base is 0.5
        score = calculate_confidence(results, "unknown_param")
        assert score == 0.5

    def test_multiple_files_base_confidence(self):
        """Test base confidence increases with more files."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
            )
            for _ in range(5)
        ]

        # 5+ files base is 0.9
        score = calculate_confidence(results, "unknown_param")
        assert score == 0.9

    def test_two_files_base_confidence(self):
        """Test base confidence with 2 files."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
            )
            for _ in range(2)
        ]

        score = calculate_confidence(results, "unknown_param")
        assert score == 0.7


class TestInstrumentConfidence:
    """Tests for _instrument_confidence function."""

    def test_consistent_instrument(self):
        """Test confidence when all files have same instrument."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive HF"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
            )
            for _ in range(3)
        ]

        base = 0.8
        score = _instrument_confidence(results, base)

        # Should be base + 0.1 (max 1.0)
        assert score == 0.9

    def test_all_unknown_instrument(self):
        """Test confidence when all instruments are unknown."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "unknown"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
            )
            for _ in range(3)
        ]

        base = 0.8
        score = _instrument_confidence(results, base)

        # All unknown: base * 0.8
        assert score == pytest.approx(0.64)

    def test_inconsistent_instruments(self):
        """Test confidence when instruments differ."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
            ),
            MockAnalysisResult(
                instrument_model={"name": "Orbitrap Fusion"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
            ),
        ]

        base = 0.7
        score = _instrument_confidence(results, base)

        # Inconsistent: base * 0.6
        assert score == pytest.approx(0.42)


class TestFragmentationConfidence:
    """Tests for _fragmentation_confidence function."""

    def test_consistent_fragmentation(self):
        """Test confidence when all files have same fragmentation."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
            )
            for _ in range(3)
        ]

        base = 0.8
        score = _fragmentation_confidence(results, base)

        assert score == 0.9

    def test_multiple_fragmentation_types(self):
        """Test confidence when multiple fragmentation types detected across files."""
        # Need multiple results with "multiple" type to trigger the penalty
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="multiple",
                precursor_tolerance_ppm=10.0,
            ),
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
            ),
        ]

        base = 0.7
        score = _fragmentation_confidence(results, base)

        # Multiple types across files: base * 0.5
        assert score == pytest.approx(0.35)

    def test_unknown_fragmentation(self):
        """Test confidence when fragmentation is unknown."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="unknown",
                precursor_tolerance_ppm=10.0,
            ),
        ]

        base = 0.5
        score = _fragmentation_confidence(results, base)

        # Unknown: base * 0.7
        assert score == pytest.approx(0.35)


class TestToleranceConfidence:
    """Tests for _tolerance_confidence function."""

    def test_no_values_detected(self):
        """Test confidence when no tolerance values detected."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=None,
            ),
        ]

        base = 0.5
        score = _tolerance_confidence(results, "precursor", base)

        assert score == 0.3

    def test_single_value(self):
        """Test confidence with single tolerance value."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
            ),
        ]

        base = 0.5
        score = _tolerance_confidence(results, "precursor", base)

        # Single value: base * 0.7
        assert score == pytest.approx(0.35)

    def test_consistent_tolerance_values(self):
        """Test confidence when tolerance values are consistent."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
            ),
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.2,
            ),
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=9.8,
            ),
        ]

        base = 0.8
        score = _tolerance_confidence(results, "precursor", base)

        # Very consistent (CV < 0.1): base + 0.1
        assert score == pytest.approx(0.9)

    def test_variable_tolerance_values(self):
        """Test confidence when tolerance values vary significantly."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=5.0,
            ),
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=20.0,
            ),
        ]

        base = 0.7
        score = _tolerance_confidence(results, "precursor", base)

        # High variability: base * 0.6
        assert score == pytest.approx(0.42)

    def test_fragment_tolerance_ppm(self):
        """Test fragment tolerance with ppm values."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
                fragment_tolerance_ppm=15.0,
            ),
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
                fragment_tolerance_ppm=14.5,
            ),
        ]

        base = 0.7
        score = _tolerance_confidence(results, "fragment", base)

        assert score > 0.5

    def test_fragment_tolerance_da(self):
        """Test fragment tolerance with Da values."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
                fragment_tolerance_da=0.02,
            ),
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
                fragment_tolerance_da=0.02,
            ),
        ]

        base = 0.7
        score = _tolerance_confidence(results, "fragment", base)

        assert score > 0.5


class TestScoreToLabel:
    """Tests for _score_to_label function."""

    def test_high_score(self):
        """Test HIGH label for scores >= 0.9."""
        assert _score_to_label(0.9) == "HIGH"
        assert _score_to_label(1.0) == "HIGH"

    def test_medium_score(self):
        """Test MEDIUM label for scores >= 0.7."""
        assert _score_to_label(0.7) == "MEDIUM"
        assert _score_to_label(0.85) == "MEDIUM"

    def test_low_score(self):
        """Test LOW label for scores >= 0.5."""
        assert _score_to_label(0.5) == "LOW"
        assert _score_to_label(0.65) == "LOW"

    def test_very_low_score(self):
        """Test VERY LOW label for scores < 0.5."""
        assert _score_to_label(0.4) == "VERY LOW"
        assert _score_to_label(0.0) == "VERY LOW"


class TestGetConfidenceFactors:
    """Tests for get_confidence_factors function."""

    def test_returns_all_parameters(self):
        """Test that confidence factors are returned for all parameters."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
                fragment_tolerance_ppm=15.0,
            ),
        ]

        factors = get_confidence_factors(results)

        assert "instrument" in factors
        assert "fragmentation" in factors
        assert "precursor_tolerance" in factors
        assert "fragment_tolerance" in factors

    def test_factor_structure(self):
        """Test structure of returned confidence factors."""
        results = [
            MockAnalysisResult(
                instrument_model={"name": "Q Exactive"},
                fragmentation_type="HR_HCD",
                precursor_tolerance_ppm=10.0,
            ),
        ]

        factors = get_confidence_factors(results)

        for param, factor in factors.items():
            assert "score" in factor
            assert "label" in factor
            assert "files_analyzed" in factor
            assert isinstance(factor["score"], float)
            assert factor["label"] in ["HIGH", "MEDIUM", "LOW", "VERY LOW"]
            assert factor["files_analyzed"] == 1
