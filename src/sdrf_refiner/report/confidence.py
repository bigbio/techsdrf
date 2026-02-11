"""
Confidence scoring utilities for SDRF refinement.
"""

from typing import Dict, List

from sdrf_refiner.analyzer.ms_analyzer import AnalysisResult


def calculate_confidence(
    analysis_results: List[AnalysisResult],
    parameter_name: str,
) -> float:
    """
    Calculate confidence score for a parameter based on analysis results.

    Factors considered:
    - Number of files analyzed
    - Consistency across files
    - Instrument category reliability
    - Fragmentation type clarity

    Args:
        analysis_results: List of per-file analysis results
        parameter_name: Name of parameter to score

    Returns:
        Confidence score between 0.0 and 1.0
    """
    if not analysis_results:
        return 0.0

    num_files = len(analysis_results)

    # Base confidence from number of files
    if num_files >= 5:
        base_confidence = 0.9
    elif num_files >= 3:
        base_confidence = 0.8
    elif num_files >= 2:
        base_confidence = 0.7
    else:
        base_confidence = 0.5

    # Adjust based on parameter type
    if parameter_name == "instrument":
        return _instrument_confidence(analysis_results, base_confidence)
    elif parameter_name == "fragmentation":
        return _fragmentation_confidence(analysis_results, base_confidence)
    elif parameter_name == "precursor_tolerance":
        return _tolerance_confidence(analysis_results, "precursor", base_confidence)
    elif parameter_name == "fragment_tolerance":
        return _tolerance_confidence(analysis_results, "fragment", base_confidence)
    else:
        return base_confidence


def _instrument_confidence(
    results: List[AnalysisResult], base: float
) -> float:
    """Calculate confidence for instrument detection."""
    # Count unique instruments
    instruments = [r.instrument_model.get("name", "unknown") for r in results]
    unique = set(instruments)

    if len(unique) == 1 and "unknown" not in unique:
        return min(base + 0.1, 1.0)  # All agree, high confidence
    elif len(unique) == 1:
        return base * 0.8  # All unknown
    else:
        return base * 0.6  # Inconsistent


def _fragmentation_confidence(
    results: List[AnalysisResult], base: float
) -> float:
    """Calculate confidence for fragmentation type detection."""
    frag_types = [r.fragmentation_type for r in results]
    unique = set(frag_types)

    if len(unique) == 1 and "unknown" not in unique:
        return min(base + 0.1, 1.0)
    elif "multiple" in unique:
        return base * 0.5  # Multiple types detected
    elif len(unique) == 1:
        return base * 0.7
    else:
        return base * 0.6


def _tolerance_confidence(
    results: List[AnalysisResult],
    tolerance_type: str,
    base: float,
) -> float:
    """Calculate confidence for tolerance values."""
    import numpy as np

    if tolerance_type == "precursor":
        values = [
            r.precursor_tolerance_ppm
            for r in results
            if r.precursor_tolerance_ppm is not None
        ]
    else:
        values = [
            r.fragment_tolerance_ppm or r.fragment_tolerance_da
            for r in results
            if r.fragment_tolerance_ppm is not None or r.fragment_tolerance_da is not None
        ]

    if not values:
        return 0.3  # No values detected

    if len(values) == 1:
        return base * 0.7

    # Check consistency (coefficient of variation)
    mean_val = np.mean(values)
    std_val = np.std(values)

    if mean_val == 0:
        return base * 0.5

    cv = std_val / mean_val

    if cv < 0.1:
        return min(base + 0.1, 1.0)  # Very consistent
    elif cv < 0.25:
        return base
    elif cv < 0.5:
        return base * 0.8
    else:
        return base * 0.6  # High variability


def get_confidence_factors(results: List[AnalysisResult]) -> Dict[str, Dict]:
    """
    Get detailed confidence factors for reporting.

    Args:
        results: List of analysis results

    Returns:
        Dictionary with confidence breakdown by parameter
    """
    factors = {}

    for param in ["instrument", "fragmentation", "precursor_tolerance", "fragment_tolerance"]:
        score = calculate_confidence(results, param)
        factors[param] = {
            "score": score,
            "label": _score_to_label(score),
            "files_analyzed": len(results),
        }

    return factors


def _score_to_label(score: float) -> str:
    """Convert score to label."""
    if score >= 0.9:
        return "HIGH"
    elif score >= 0.7:
        return "MEDIUM"
    elif score >= 0.5:
        return "LOW"
    else:
        return "VERY LOW"
