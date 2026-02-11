"""
Tests for ParameterComparer - comparison logic between SDRF and detected values.
"""

import pytest

from sdrf_refiner.sdrf.comparer import (
    ParameterComparer,
    ParameterComparison,
    ComparisonStatus,
)


class TestCompareInstrument:
    """Tests for instrument comparison."""

    def test_match_exact(self):
        """Test exact instrument match."""
        comparer = ParameterComparer()
        result = comparer.compare_instrument(
            {"name": "Q Exactive HF", "accession": "MS:1002523"},
            {"name": "Q Exactive HF", "accession": "MS:1002523"},
        )

        assert result.status == ComparisonStatus.MATCH
        assert result.recommendation is None

    def test_match_partial_name(self):
        """Test partial name match (fuzzy matching)."""
        comparer = ParameterComparer()
        result = comparer.compare_instrument(
            {"name": "Q Exactive", "accession": ""},
            {"name": "Q Exactive HF", "accession": "MS:1002523"},
        )

        assert result.status == ComparisonStatus.MATCH

    def test_mismatch(self):
        """Test instrument mismatch."""
        comparer = ParameterComparer()
        result = comparer.compare_instrument(
            {"name": "Orbitrap Fusion", "accession": "MS:1002416"},
            {"name": "Q Exactive HF", "accession": "MS:1002523"},
        )

        assert result.status == ComparisonStatus.MISMATCH
        assert "Q Exactive HF" in result.recommendation

    def test_missing_sdrf(self):
        """Test when SDRF has no instrument specified."""
        comparer = ParameterComparer()
        result = comparer.compare_instrument(
            {"name": "", "accession": ""},
            {"name": "Q Exactive HF", "accession": "MS:1002523"},
        )

        assert result.status == ComparisonStatus.MISSING_SDRF
        assert "Add instrument" in result.recommendation

    def test_missing_detected(self):
        """Test when detection returns unknown."""
        comparer = ParameterComparer()
        result = comparer.compare_instrument(
            {"name": "Q Exactive HF", "accession": "MS:1002523"},
            {"name": "unknown"},
        )

        assert result.status == ComparisonStatus.MISSING_DETECTED

    def test_detected_as_string(self):
        """Test when detected value is a plain string."""
        comparer = ParameterComparer()
        result = comparer.compare_instrument(
            {"name": "Q Exactive HF", "accession": "MS:1002523"},
            "Q Exactive HF",
        )

        assert result.status == ComparisonStatus.MATCH

    def test_none_sdrf_value(self):
        """Test when SDRF value is None."""
        comparer = ParameterComparer()
        result = comparer.compare_instrument(
            None,
            {"name": "Q Exactive HF", "accession": "MS:1002523"},
        )

        assert result.status == ComparisonStatus.MISSING_SDRF


class TestCompareFragmentation:
    """Tests for fragmentation/dissociation comparison."""

    def test_match_hcd(self):
        """Test HCD fragmentation match."""
        comparer = ParameterComparer()
        result = comparer.compare_fragmentation(
            {"name": "HCD", "accession": "MS:1000422"},
            "HR_HCD",
        )

        assert result.status == ComparisonStatus.MATCH

    def test_mismatch_cid_vs_hcd(self):
        """Test mismatch between CID and HCD."""
        comparer = ParameterComparer()
        result = comparer.compare_fragmentation(
            {"name": "CID", "accession": "MS:1000133"},
            "HR_HCD",
        )

        assert result.status == ComparisonStatus.MISMATCH
        assert "HCD" in result.recommendation

    def test_missing_sdrf(self):
        """Test when SDRF has no dissociation method."""
        comparer = ParameterComparer()
        result = comparer.compare_fragmentation(
            {"name": "", "accession": ""},
            "HR_HCD",
        )

        assert result.status == ComparisonStatus.MISSING_SDRF
        assert "Add dissociation method" in result.recommendation

    def test_unknown_detected_type(self):
        """Test with unknown detected fragmentation type."""
        comparer = ParameterComparer()
        result = comparer.compare_fragmentation(
            {"name": "HCD", "accession": "MS:1000422"},
            "UNKNOWN_TYPE",
        )

        assert result.status == ComparisonStatus.MISSING_DETECTED

    def test_ms3_dissociation_param_name(self):
        """Test MS3 dissociation with custom param name."""
        comparer = ParameterComparer()
        result = comparer.compare_fragmentation(
            {"name": "", "accession": ""},
            "HR_HCD",
            param_name="ms3_dissociation_method",
        )

        assert result.status == ComparisonStatus.MISSING_SDRF
        assert "MS3 dissociation method" in result.recommendation


class TestCompareTolerance:
    """Tests for tolerance comparison."""

    def test_match_within_threshold(self):
        """Test tolerance match within threshold (20%)."""
        comparer = ParameterComparer()
        result = comparer.compare_tolerance(
            "precursor_mass_tolerance",
            {"value": 10.0, "unit": "ppm"},
            11.0,
            "ppm",
        )

        assert result.status == ComparisonStatus.MATCH

    def test_mismatch_above_threshold(self):
        """Test tolerance mismatch when difference exceeds threshold."""
        comparer = ParameterComparer()
        result = comparer.compare_tolerance(
            "precursor_mass_tolerance",
            {"value": 10.0, "unit": "ppm"},
            25.0,
            "ppm",
        )

        assert result.status == ComparisonStatus.MISMATCH
        assert "25" in result.recommendation

    def test_improved_more_stringent(self):
        """Test when detected value is more stringent (lower)."""
        comparer = ParameterComparer()
        result = comparer.compare_tolerance(
            "precursor_mass_tolerance",
            {"value": 20.0, "unit": "ppm"},
            8.0,
            "ppm",
        )

        assert result.status == ComparisonStatus.IMPROVED
        assert "8" in result.recommendation

    def test_missing_sdrf(self):
        """Test when SDRF has no tolerance specified."""
        comparer = ParameterComparer()
        result = comparer.compare_tolerance(
            "precursor_mass_tolerance",
            {"value": None, "unit": ""},
            10.0,
            "ppm",
        )

        assert result.status == ComparisonStatus.MISSING_SDRF
        assert "Add precursor_mass_tolerance" in result.recommendation

    def test_missing_detected(self):
        """Test when detection returns None."""
        comparer = ParameterComparer()
        result = comparer.compare_tolerance(
            "fragment_mass_tolerance",
            {"value": 10.0, "unit": "ppm"},
            None,
            "ppm",
        )

        assert result.status == ComparisonStatus.MISSING_DETECTED

    def test_unit_mismatch_cross_unit_comparison(self):
        """Test cross-unit comparison: SDRF=10ppm vs detected=0.02 Da (~40ppm at 500 m/z).

        The comparer now converts both to ppm for comparison.
        0.02 Da at 500 m/z = 40 ppm, which is wider than 10 ppm → MISMATCH.
        """
        comparer = ParameterComparer(reference_mz=500.0)
        result = comparer.compare_tolerance(
            "fragment_mass_tolerance",
            {"value": 10.0, "unit": "ppm"},
            0.02,
            "Da",
        )

        assert result.status == ComparisonStatus.MISMATCH
        assert result.details.get("tolerance_source") == "empirical"
        # The comparison was done in ppm: 10 ppm vs 40 ppm
        assert result.details.get("comparison_ppm", {}).get("sdrf") == 10.0
        assert result.details.get("comparison_ppm", {}).get("detected") == 40.0

    def test_hard_divergence_always_flags(self):
        """Test >10x divergence always flagged, even when methods agree.

        SDRF=0.5 Da (~1000 ppm at 500 m/z) vs detected=20 ppm → 50x divergence.
        Exceeds hard limit (10x), so always flagged.
        """
        comparer = ParameterComparer(reference_mz=500.0)
        # Without methods_agree
        result = comparer.compare_tolerance(
            "fragment_mass_tolerance",
            {"value": 0.5, "unit": "Da"},
            20.0,
            "ppm",
        )
        assert result.status == ComparisonStatus.MISSING_DETECTED
        assert result.details.get("divergence_factor") == 50.0

        # Even with methods_agree=True, still flagged (>10x hard limit)
        result = comparer.compare_tolerance(
            "fragment_mass_tolerance",
            {"value": 0.5, "unit": "Da"},
            20.0,
            "ppm",
            methods_agree=True,
        )
        assert result.status == ComparisonStatus.MISSING_DETECTED

    def test_soft_divergence_single_method_flags(self):
        """Test 5x-10x divergence flagged when methods don't agree.

        SDRF=0.04 Da (~80 ppm at 500 m/z) vs detected=10 ppm → 8x divergence.
        Between soft (5x) and hard (10x) limits, single method → flagged.
        """
        comparer = ParameterComparer(reference_mz=500.0)
        result = comparer.compare_tolerance(
            "fragment_mass_tolerance",
            {"value": 0.04, "unit": "Da"},
            10.0,
            "ppm",
            methods_agree=False,
        )
        assert result.status == ComparisonStatus.MISSING_DETECTED

    def test_soft_divergence_methods_agree_overrides(self):
        """Test 5x-10x divergence allowed when methods agree.

        SDRF=0.04 Da (~80 ppm at 500 m/z) vs detected=10 ppm → 8x divergence.
        Between soft (5x) and hard (10x) limits, methods agree → override.
        """
        comparer = ParameterComparer(reference_mz=500.0)
        result = comparer.compare_tolerance(
            "fragment_mass_tolerance",
            {"value": 0.04, "unit": "Da"},
            10.0,
            "ppm",
            methods_agree=True,
        )
        assert result.status == ComparisonStatus.IMPROVED

    def test_moderate_divergence_improved(self):
        """Test <5x divergence → IMPROVED regardless.

        SDRF=0.05 Da (~100 ppm at 500 m/z) vs detected=40 ppm → 2.5x divergence.
        Within soft limit, so IMPROVED even without methods_agree.
        """
        comparer = ParameterComparer(reference_mz=500.0)
        result = comparer.compare_tolerance(
            "fragment_mass_tolerance",
            {"value": 0.05, "unit": "Da"},
            40.0,
            "ppm",
        )
        assert result.status == ComparisonStatus.IMPROVED

    def test_heuristic_only_keeps_sdrf(self):
        """Test that heuristic-only tolerance does NOT change SDRF."""
        comparer = ParameterComparer(reference_mz=500.0)
        result = comparer.compare_tolerance(
            "fragment_mass_tolerance",
            {"value": 0.2, "unit": "Da"},
            20.0,
            "ppm",
            tolerance_source="heuristic",
            instrument_name="Q Exactive",
        )

        # Heuristic-only should NOT trigger MISMATCH/IMPROVED
        assert result.status == ComparisonStatus.MISSING_DETECTED
        assert result.details.get("heuristic_only") is True
        assert "recommend" in result.recommendation.lower()


class TestCompareCharge:
    """Tests for charge state comparison."""

    def test_match(self):
        """Test charge match."""
        comparer = ParameterComparer()
        result = comparer.compare_charge("ms_min_charge", 2, 2)

        assert result.status == ComparisonStatus.MATCH

    def test_mismatch(self):
        """Test charge mismatch."""
        comparer = ParameterComparer()
        result = comparer.compare_charge("ms_max_charge", 6, 8)

        assert result.status == ComparisonStatus.MISMATCH
        assert "8" in result.recommendation

    def test_missing_sdrf(self):
        """Test when SDRF has no charge specified."""
        comparer = ParameterComparer()
        result = comparer.compare_charge("ms_min_charge", None, 2)

        assert result.status == ComparisonStatus.MISSING_SDRF

    def test_missing_detected(self):
        """Test when detection returns None."""
        comparer = ParameterComparer()
        result = comparer.compare_charge("ms_max_charge", 6, None)

        assert result.status == ComparisonStatus.MISSING_DETECTED


class TestCompareAll:
    """Tests for compare_all method."""

    def test_compare_all_basic(self):
        """Test comparing all parameters."""
        comparer = ParameterComparer()

        sdrf_params = {
            "instrument": {"name": "Q Exactive HF", "accession": "MS:1002523"},
            "dissociation": {"name": "HCD", "accession": "MS:1000422"},
            "precursor_tolerance": {"value": 10.0, "unit": "ppm"},
            "fragment_tolerance": {"value": 0.02, "unit": "Da"},
            "ms_min_charge": 2,
            "ms_max_charge": 6,
        }

        detected_params = {
            "instrument_model": {"name": "Q Exactive HF", "accession": "MS:1002523"},
            "fragmentation_type": "HR_HCD",
            "precursor_tolerance_ppm": 10.5,
            "fragment_tolerance_da": 0.02,
            "ms_min_charge": 2,
            "ms_max_charge": 6,
        }

        comparisons = comparer.compare_all(sdrf_params, detected_params)

        # Should have comparisons for instrument, fragmentation, precursor, fragment, min_charge, max_charge
        assert len(comparisons) >= 6

        # Most should match
        matches = [c for c in comparisons if c.status == ComparisonStatus.MATCH]
        assert len(matches) >= 4


class TestGetRefinements:
    """Tests for get_refinements method."""

    def test_get_refinements_filters_matches(self):
        """Test that get_refinements excludes matches."""
        comparer = ParameterComparer()
        comparer.comparisons = [
            ParameterComparison(
                parameter_name="instrument",
                sdrf_value={},
                detected_value={},
                status=ComparisonStatus.MATCH,
            ),
            ParameterComparison(
                parameter_name="tolerance",
                sdrf_value={},
                detected_value={},
                status=ComparisonStatus.MISMATCH,
                recommendation="Update tolerance",
            ),
        ]

        refinements = comparer.get_refinements()

        assert len(refinements) == 1
        assert refinements[0].parameter_name == "tolerance"

    def test_get_refinements_includes_missing_sdrf(self):
        """Test that get_refinements includes MISSING_SDRF status."""
        comparer = ParameterComparer()
        comparer.comparisons = [
            ParameterComparison(
                parameter_name="instrument",
                sdrf_value=None,
                detected_value={"name": "Q Exactive"},
                status=ComparisonStatus.MISSING_SDRF,
                recommendation="Add instrument",
            ),
        ]

        refinements = comparer.get_refinements()

        assert len(refinements) == 1
        assert refinements[0].status == ComparisonStatus.MISSING_SDRF

    def test_get_refinements_includes_improved(self):
        """Test that get_refinements includes IMPROVED status."""
        comparer = ParameterComparer()
        comparer.comparisons = [
            ParameterComparison(
                parameter_name="tolerance",
                sdrf_value={"value": 20},
                detected_value={"value": 10},
                status=ComparisonStatus.IMPROVED,
                recommendation="Update tolerance",
            ),
        ]

        refinements = comparer.get_refinements()

        assert len(refinements) == 1
        assert refinements[0].status == ComparisonStatus.IMPROVED


class TestCompareNumericMz:
    """Tests for compare_numeric_mz (scan windows, isolation width)."""

    def test_match(self):
        comparer = ParameterComparer()
        comp = comparer.compare_numeric_mz("scan_window_lower", 350.0, 350.0)
        assert comp.status == ComparisonStatus.MATCH

    def test_match_within_1_percent(self):
        comparer = ParameterComparer()
        comp = comparer.compare_numeric_mz("scan_window_upper", 1600.0, 1605.0)
        assert comp.status == ComparisonStatus.MATCH

    def test_mismatch(self):
        comparer = ParameterComparer()
        comp = comparer.compare_numeric_mz("scan_window_lower", 400.0, 350.0)
        assert comp.status == ComparisonStatus.MISMATCH

    def test_missing_sdrf(self):
        comparer = ParameterComparer()
        comp = comparer.compare_numeric_mz("scan_window_lower", None, 350.0)
        assert comp.status == ComparisonStatus.MISSING_SDRF
        assert "350" in comp.recommendation

    def test_missing_detected(self):
        comparer = ParameterComparer()
        comp = comparer.compare_numeric_mz("scan_window_lower", 350.0, None)
        assert comp.status == ComparisonStatus.MISSING_DETECTED


class TestCompareAllScanWindows:
    """Tests for scan window comparisons in compare_all."""

    def test_scan_windows_included(self):
        comparer = ParameterComparer()
        sdrf_params = {}
        detected_params = {
            "ms1_scan_range": (350.0, 1600.0),
            "isolation_window_da": 2.0,
        }
        results = comparer.compare_all(sdrf_params, detected_params)
        param_names = [c.parameter_name for c in results]
        assert "scan_window_lower" in param_names
        assert "scan_window_upper" in param_names
        assert "isolation_window_width" in param_names

    def test_scan_windows_missing_sdrf(self):
        comparer = ParameterComparer()
        sdrf_params = {}
        detected_params = {
            "ms1_scan_range": (350.0, 1600.0),
            "isolation_window_da": 25.0,
        }
        results = comparer.compare_all(sdrf_params, detected_params)
        sw_lower = next(c for c in results if c.parameter_name == "scan_window_lower")
        sw_upper = next(c for c in results if c.parameter_name == "scan_window_upper")
        iw = next(c for c in results if c.parameter_name == "isolation_window_width")

        assert sw_lower.status == ComparisonStatus.MISSING_SDRF
        assert sw_lower.detected_value == 350.0
        assert sw_upper.status == ComparisonStatus.MISSING_SDRF
        assert sw_upper.detected_value == 1600.0
        assert iw.status == ComparisonStatus.MISSING_SDRF
        assert iw.detected_value == 25.0

    def test_scan_windows_not_included_when_not_detected(self):
        comparer = ParameterComparer()
        sdrf_params = {}
        detected_params = {}
        results = comparer.compare_all(sdrf_params, detected_params)
        param_names = [c.parameter_name for c in results]
        assert "scan_window_lower" not in param_names
        assert "scan_window_upper" not in param_names
        assert "isolation_window_width" not in param_names
