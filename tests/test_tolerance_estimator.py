"""
Tests for tolerance estimation module.
"""

import pytest
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import numpy as np

from sdrf_refiner.analyzer.tolerance_estimator import (
    ppm_to_da,
    da_to_ppm,
    convert_tolerance,
    normalize_tolerance,
    ToleranceResult,
    ConsensusResult,
    GaussianEstimator,
    CruxEstimator,
    ToleranceEstimator,
)


class TestUnitConversion:
    """Tests for unit conversion utilities."""

    def test_ppm_to_da(self):
        """Test ppm to Da conversion."""
        # At m/z 500, 10 ppm = 0.005 Da
        result = ppm_to_da(10.0, 500.0)
        assert result == pytest.approx(0.005)

        # At m/z 1000, 10 ppm = 0.01 Da
        result = ppm_to_da(10.0, 1000.0)
        assert result == pytest.approx(0.01)

    def test_da_to_ppm(self):
        """Test Da to ppm conversion."""
        # At m/z 500, 0.005 Da = 10 ppm
        result = da_to_ppm(0.005, 500.0)
        assert result == pytest.approx(10.0)

        # At m/z 1000, 0.01 Da = 10 ppm
        result = da_to_ppm(0.01, 1000.0)
        assert result == pytest.approx(10.0)

    def test_da_to_ppm_zero_mz(self):
        """Test Da to ppm with zero m/z returns 0."""
        result = da_to_ppm(0.005, 0.0)
        assert result == 0.0

    def test_normalize_tolerance_ppm(self):
        """Test normalize tolerance from ppm."""
        ppm_val, da_val = normalize_tolerance(10.0, 'ppm', reference_mz=500.0)
        assert ppm_val == 10.0
        assert da_val == pytest.approx(0.005)

    def test_normalize_tolerance_da(self):
        """Test normalize tolerance from Da."""
        ppm_val, da_val = normalize_tolerance(0.02, 'Da', reference_mz=500.0)
        assert da_val == 0.02
        assert ppm_val == pytest.approx(40.0)

    def test_roundtrip_conversion(self):
        """Test that conversions roundtrip correctly."""
        original_ppm = 15.0
        mz = 750.0

        da = ppm_to_da(original_ppm, mz)
        back_to_ppm = da_to_ppm(da, mz)

        assert back_to_ppm == pytest.approx(original_ppm)

    def test_convert_tolerance_ppm_to_da(self):
        """Test generic convert_tolerance: ppm -> Da."""
        result = convert_tolerance(10.0, "ppm", "Da", reference_mz=500.0)
        assert result == pytest.approx(0.005)

    def test_convert_tolerance_da_to_ppm(self):
        """Test generic convert_tolerance: Da -> ppm."""
        result = convert_tolerance(0.02, "Da", "ppm", reference_mz=500.0)
        assert result == pytest.approx(40.0)

    def test_convert_tolerance_unit_aliases(self):
        """Test convert_tolerance accepts Th and m/z as Da aliases."""
        # Th -> ppm
        assert convert_tolerance(0.02, "Th", "ppm", reference_mz=500.0) == pytest.approx(40.0)
        # m/z -> ppm
        assert convert_tolerance(0.02, "m/z", "ppm", reference_mz=500.0) == pytest.approx(40.0)
        # ppm -> Th (Da)
        assert convert_tolerance(40.0, "ppm", "Th", reference_mz=500.0) == pytest.approx(0.02)

    def test_convert_tolerance_same_unit_identity(self):
        """Test convert_tolerance returns identity for same-unit conversion."""
        assert convert_tolerance(10.0, "ppm", "ppm") == 10.0
        assert convert_tolerance(0.02, "Da", "Da") == 0.02

    def test_convert_tolerance_unknown_unit_raises(self):
        """Test convert_tolerance raises for unknown unit."""
        with pytest.raises(ValueError, match="Unknown mass tolerance unit"):
            convert_tolerance(10.0, "unknown", "Da")
        with pytest.raises(ValueError, match="Unknown mass tolerance unit"):
            convert_tolerance(10.0, "ppm", "kg")

    def test_normalize_tolerance_with_th_alias(self):
        """Test normalize_tolerance works with Th (Da) alias."""
        ppm_val, da_val = normalize_tolerance(0.02, "Th", reference_mz=500.0)
        assert da_val == 0.02
        assert ppm_val == pytest.approx(40.0)


class TestToleranceResult:
    """Tests for ToleranceResult dataclass."""

    def test_has_precursor(self):
        """Test has_precursor method."""
        result = ToleranceResult(method="test", precursor_ppm=10.0)
        assert result.has_precursor() is True

        result = ToleranceResult(method="test", precursor_da=0.01)
        assert result.has_precursor() is True

        result = ToleranceResult(method="test")
        assert result.has_precursor() is False

    def test_has_fragment(self):
        """Test has_fragment method."""
        result = ToleranceResult(method="test", fragment_ppm=20.0)
        assert result.has_fragment() is True

        result = ToleranceResult(method="test", fragment_da=0.02)
        assert result.has_fragment() is True

        result = ToleranceResult(method="test")
        assert result.has_fragment() is False

    def test_get_precursor_ppm_direct(self):
        """Test get_precursor_ppm when ppm is set directly."""
        result = ToleranceResult(method="test", precursor_ppm=10.0)
        assert result.get_precursor_ppm() == 10.0

    def test_get_precursor_ppm_from_da(self):
        """Test get_precursor_ppm conversion from Da."""
        result = ToleranceResult(method="test", precursor_da=0.005)
        ppm = result.get_precursor_ppm(reference_mz=500.0)
        assert ppm == pytest.approx(10.0)

    def test_get_fragment_ppm_direct(self):
        """Test get_fragment_ppm when ppm is set directly."""
        result = ToleranceResult(method="test", fragment_ppm=20.0)
        assert result.get_fragment_ppm() == 20.0

    def test_get_fragment_ppm_from_da(self):
        """Test get_fragment_ppm conversion from Da."""
        result = ToleranceResult(method="test", fragment_da=0.01)
        ppm = result.get_fragment_ppm(reference_mz=500.0)
        assert ppm == pytest.approx(20.0)

    def test_get_precursor_da_direct(self):
        """Test get_precursor_da when Da is set directly."""
        result = ToleranceResult(method="test", precursor_da=0.005)
        assert result.get_precursor_da() == 0.005

    def test_get_precursor_da_from_ppm(self):
        """Test get_precursor_da conversion from ppm."""
        result = ToleranceResult(method="test", precursor_ppm=10.0)
        da = result.get_precursor_da(reference_mz=500.0)
        assert da == pytest.approx(0.005)

    def test_get_fragment_da_direct(self):
        """Test get_fragment_da when Da is set directly."""
        result = ToleranceResult(method="test", fragment_da=0.02)
        assert result.get_fragment_da() == 0.02

    def test_get_fragment_da_from_ppm(self):
        """Test get_fragment_da conversion from ppm."""
        result = ToleranceResult(method="test", fragment_ppm=40.0)
        da = result.get_fragment_da(reference_mz=500.0)
        assert da == pytest.approx(0.02)


class TestConsensusResult:
    """Tests for ConsensusResult dataclass."""

    def test_format_report_basic(self):
        """Test basic report formatting."""
        consensus = ConsensusResult()
        consensus.gaussian = ToleranceResult(method="gaussian", precursor_ppm=10.0)
        consensus.recommended_precursor_ppm = 10.0
        consensus.recommended_precursor_da = 0.005
        consensus.precursor_confidence = 0.9

        report = consensus.format_report()

        assert "TOLERANCE ESTIMATION REPORT" in report
        assert "PRECURSOR MASS TOLERANCE" in report
        assert "Gaussian:" in report
        assert "RECOMMENDED:" in report

    def test_format_report_all_methods(self):
        """Test report with Gaussian and Crux methods."""
        consensus = ConsensusResult()
        consensus.gaussian = ToleranceResult(
            method="gaussian",
            precursor_ppm=10.0,
            confidence=0.9,
            details={'precursor_fit_quality': 'good'}
        )
        consensus.crux = ToleranceResult(
            method="crux",
            precursor_ppm=12.0,
            fragment_da=0.02,
            confidence=0.85
        )
        consensus.recommended_precursor_ppm = 11.0
        consensus.recommended_precursor_da = 0.0055
        consensus.recommended_fragment_ppm = 40.0  # 0.02 Da at m/z 500 ≈ 40 ppm
        consensus.recommended_fragment_da = 0.02
        consensus.precursor_confidence = 0.85
        consensus.fragment_confidence = 0.82

        report = consensus.format_report()

        assert "Gaussian:" in report
        assert "Crux:" in report
        assert "good" in report  # fit quality
        assert "0.02" in report or "0.0200" in report  # fragment Da


class TestGaussianEstimator:
    """Tests for GaussianEstimator class (RunAssessor-backed)."""

    def test_init(self):
        """Test initialization."""
        estimator = GaussianEstimator()
        assert estimator is not None

    @patch('sdrf_refiner.analyzer.tolerance_estimator.GaussianEstimator.estimate')
    def test_estimate_returns_tolerance_result(self, mock_estimate):
        """Test that estimate returns a valid ToleranceResult."""
        mock_estimate.return_value = ToleranceResult(
            method="gaussian",
            precursor_ppm=10.0,
            confidence=0.8,
            details={"source": "runassessor"},
        )
        estimator = GaussianEstimator()
        result = estimator.estimate(Path("/fake/file.mzML"))
        assert result.method == "gaussian"
        assert result.precursor_ppm == 10.0


class TestCruxEstimator:
    """Tests for CruxEstimator class."""

    def test_init(self):
        """Test initialization."""
        estimator = CruxEstimator(charges="2,3", min_peak_pairs=30)
        assert estimator.charges == "2,3"
        assert estimator.min_peak_pairs == 30

    @patch('subprocess.run')
    def test_is_available_true(self, mock_run):
        """Test is_available when crux is installed."""
        mock_run.return_value = Mock(returncode=0)

        estimator = CruxEstimator()
        assert estimator.is_available() is True

    @patch('subprocess.run')
    def test_is_available_false(self, mock_run):
        """Test is_available when crux is not installed."""
        mock_run.side_effect = FileNotFoundError()

        estimator = CruxEstimator()
        assert estimator.is_available() is False

    def test_parse_crux_output(self):
        """Test parsing Crux param-medic output."""
        estimator = CruxEstimator()
        result = ToleranceResult(method="crux")

        # Mock log content
        import tempfile
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write("INFO: Precursor error: 15.50 ppm\n")
            f.write("INFO: Fragment bin size: 0.0234 Th\n")
            log_path = Path(f.name)

        try:
            estimator._parse_crux_output(log_path, result)

            assert result.precursor_ppm == pytest.approx(15.50)
            assert result.fragment_da == pytest.approx(0.0234)
        finally:
            log_path.unlink()


class TestToleranceEstimator:
    """Tests for main ToleranceEstimator orchestrator."""

    def test_init_default(self):
        """Test default initialization."""
        estimator = ToleranceEstimator()

        assert estimator.use_gaussian is True
        assert estimator.use_crux is True
        assert estimator.gaussian_estimator is not None
        assert estimator.crux_estimator is not None

    def test_init_custom(self):
        """Test custom initialization."""
        estimator = ToleranceEstimator(
            use_gaussian=False,
            use_crux=True,
        )

        assert estimator.gaussian_estimator is None
        assert estimator.crux_estimator is not None

    def test_select_best_estimate_single(self):
        """Test best estimate selection with single method."""
        estimator = ToleranceEstimator()
        consensus = ConsensusResult()

        estimates = [(10.0, 0.9, 'gaussian')]
        value, conf, basis = estimator._select_best_estimate(estimates, consensus)

        assert value == 10.0
        assert conf == 0.9
        assert "single method" in basis

    def test_select_best_estimate_all_agree(self):
        """Test best estimate when all methods agree."""
        estimator = ToleranceEstimator()
        consensus = ConsensusResult()

        # All within 20% of each other
        estimates = [
            (10.0, 0.9, 'gaussian'),
            (11.0, 0.85, 'crux'),
        ]
        value, conf, basis = estimator._select_best_estimate(estimates, consensus)

        # Should be weighted average
        assert value == pytest.approx(10.5, rel=0.1)
        assert "weighted average" in basis
        assert "all agree" in basis

    def test_select_best_estimate_two_agree(self):
        """Test best estimate when two methods agree, one differs."""
        estimator = ToleranceEstimator()
        consensus = ConsensusResult()

        # Two agree, one outlier (simulating third method)
        estimates = [
            (10.0, 0.9, 'gaussian'),
            (11.0, 0.85, 'crux'),
            (50.0, 0.5, 'other'),  # Outlier
        ]
        value, conf, basis = estimator._select_best_estimate(estimates, consensus)

        # Should use median (robust to outlier)
        assert value == pytest.approx(11.0)
        assert "median" in basis

    def test_select_best_estimate_all_disagree(self):
        """Test best estimate when all methods disagree."""
        estimator = ToleranceEstimator()
        consensus = ConsensusResult()

        # All significantly different
        estimates = [
            (10.0, 0.9, 'gaussian'),
            (50.0, 0.85, 'crux'),
        ]
        value, conf, basis = estimator._select_best_estimate(estimates, consensus)

        # Should prefer highest confidence
        assert value == 10.0  # gaussian has highest confidence
        assert "highest confidence" in basis
        assert "gaussian" in basis

    @patch.object(GaussianEstimator, 'estimate')
    @patch.object(CruxEstimator, 'estimate')
    def test_estimate_integration(self, mock_crux, mock_gaussian):
        """Test full estimation workflow and fragment tolerance in Da."""
        mock_gaussian.return_value = ToleranceResult(
            method="gaussian",
            precursor_ppm=10.0,
            confidence=0.9,
        )
        mock_crux.return_value = ToleranceResult(
            method="crux",
            precursor_ppm=11.0,
            fragment_da=0.02,  # Crux reports fragment tolerance in Da (Th)
            confidence=0.85,
        )

        estimator = ToleranceEstimator(use_gaussian=True, use_crux=True)
        consensus = estimator.estimate([Path("/fake/file.mzML")])

        assert consensus.gaussian is not None
        assert consensus.crux is not None
        assert consensus.recommended_precursor_ppm is not None
        assert consensus.recommended_fragment_ppm is not None
        # Fragment tolerance in Da (from Crux when available)
        assert consensus.recommended_fragment_da is not None
        assert consensus.recommended_fragment_da == pytest.approx(0.02)

    def test_compute_consensus_empty(self):
        """Test consensus computation with no results."""
        estimator = ToleranceEstimator()
        consensus = ConsensusResult()

        estimator._compute_consensus(consensus)

        assert consensus.recommended_precursor_ppm is None
        assert consensus.recommended_fragment_ppm is None
