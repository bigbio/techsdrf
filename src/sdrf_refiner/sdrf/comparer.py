"""
ParameterComparer - Compare detected parameters with SDRF values.
"""

import logging
import re
from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Dict, List, Optional

from sdrf_refiner.analyzer.tolerance_estimator import convert_tolerance
from sdrf_refiner.config import (
    TOLERANCE_MATCH_THRESHOLD,
    TOLERANCE_DIVERGENCE_SOFT,
    TOLERANCE_DIVERGENCE_HARD,
    FRAGMENTATION_ONTOLOGY,
)

logger = logging.getLogger(__name__)


class ComparisonStatus(Enum):
    """Status of parameter comparison."""

    MATCH = "match"
    MISMATCH = "mismatch"
    MISSING_SDRF = "missing_in_sdrf"
    MISSING_DETECTED = "missing_detected"
    IMPROVED = "can_improve"


@dataclass
class ParameterComparison:
    """Comparison result for a single parameter."""

    parameter_name: str
    sdrf_value: Any
    detected_value: Any
    status: ComparisonStatus
    confidence: float = 1.0
    recommendation: Optional[str] = None
    details: Dict[str, Any] = field(default_factory=dict)


class ParameterComparer:
    """
    Compares empirically detected MS parameters with values in SDRF file.
    Generates recommendations for refinement.
    """

    def __init__(
        self,
        tolerance_unit: Optional[str] = None,
        reference_mz: float = 500.0,
    ):
        """
        Args:
            tolerance_unit: User-preferred unit for tolerance output ('ppm' or 'Da').
                            None = keep the unit from the detection method.
            reference_mz: Reference m/z for ppm<->Da conversions.
        """
        self.comparisons: List[ParameterComparison] = []
        self.tolerance_unit = tolerance_unit
        self.reference_mz = reference_mz

    def compare_instrument(
        self, sdrf_value: Optional[Dict], detected_value
    ) -> ParameterComparison:
        """Compare instrument values."""
        sdrf_name = (sdrf_value.get("name", "") or "").lower() if sdrf_value else ""

        # Handle detected_value being either a dict or a string
        if isinstance(detected_value, dict):
            detected_name = (detected_value.get("name", "") or "").lower()
        elif isinstance(detected_value, str):
            detected_name = detected_value.lower()
            detected_value = {"name": detected_value}  # Convert to dict for consistency
        else:
            detected_name = ""

        # Treat "unknown" as if nothing was detected - don't overwrite valid SDRF values
        if not detected_name or detected_name == "unknown":
            return ParameterComparison(
                parameter_name="instrument",
                sdrf_value=sdrf_value,
                detected_value=detected_value,
                status=ComparisonStatus.MISSING_DETECTED,
            )

        if not sdrf_name:
            return ParameterComparison(
                parameter_name="instrument",
                sdrf_value=sdrf_value,
                detected_value=detected_value,
                status=ComparisonStatus.MISSING_SDRF,
                recommendation=f"Add instrument: NT={detected_value.get('name')};AC={detected_value.get('accession', '')}",
            )

        # Check if names match (fuzzy matching for instrument names)
        if sdrf_name in detected_name or detected_name in sdrf_name:
            status = ComparisonStatus.MATCH
            recommendation = None
        else:
            status = ComparisonStatus.MISMATCH
            recommendation = f"Update instrument to: NT={detected_value.get('name')};AC={detected_value.get('accession', '')}"

        return ParameterComparison(
            parameter_name="instrument",
            sdrf_value=sdrf_value,
            detected_value=detected_value,
            status=status,
            recommendation=recommendation,
        )

    def compare_fragmentation(
        self, sdrf_value: Optional[Dict], detected_type: str,
        param_name: str = "dissociation_method"
    ) -> ParameterComparison:
        """
        Compare fragmentation/dissociation method.

        Args:
            sdrf_value: SDRF value for dissociation method
            detected_type: Detected fragmentation type (e.g., 'HR_HCD')
            param_name: Parameter name (e.g., 'dissociation_method' or 'ms3_dissociation_method')
        """
        sdrf_name = (sdrf_value.get("name", "") or "").upper() if sdrf_value else ""

        # Map detected fragmentation type to ontology term
        ont_mapping = FRAGMENTATION_ONTOLOGY.get(detected_type, {})
        detected_name = ont_mapping.get("name", "").upper()
        detected_accession = ont_mapping.get("accession", "")

        # Human-readable param name for messages
        display_name = "MS3 dissociation method" if "ms3" in param_name else "dissociation method"

        if not detected_name:
            return ParameterComparison(
                parameter_name=param_name,
                sdrf_value=sdrf_value,
                detected_value={"type": detected_type},
                status=ComparisonStatus.MISSING_DETECTED,
            )

        if not sdrf_name:
            return ParameterComparison(
                parameter_name=param_name,
                sdrf_value=sdrf_value,
                detected_value={"name": detected_name, "accession": detected_accession},
                status=ComparisonStatus.MISSING_SDRF,
                recommendation=f"Add {display_name}: NT={detected_name};AC={detected_accession}",
            )

        if sdrf_name == detected_name:
            status = ComparisonStatus.MATCH
            recommendation = None
        else:
            status = ComparisonStatus.MISMATCH
            recommendation = (
                f"Update {display_name} to: NT={detected_name};AC={detected_accession}"
            )

        return ParameterComparison(
            parameter_name=param_name,
            sdrf_value=sdrf_value,
            detected_value={"name": detected_name, "accession": detected_accession},
            status=status,
            recommendation=recommendation,
        )

    def compare_tolerance(
        self,
        parameter_name: str,
        sdrf_value: Optional[Dict],
        detected_value: Optional[float],
        detected_unit: str,
        tolerance_source: str = "empirical",
        instrument_name: str = "",
        preferred_unit: Optional[str] = None,
        methods_agree: bool = False,
    ) -> ParameterComparison:
        """
        Compare tolerance values (precursor or fragment).

        When ``tolerance_source`` is ``"heuristic"`` the detected value came
        from a generic instrument-category lookup, **not** from empirical
        estimation (Crux/Gaussian).  In that case we never change the SDRF
        value; instead we emit a review recommendation.

        ``preferred_unit`` overrides ``self.tolerance_unit`` for the output
        unit.  Pass ``"ppm"`` for precursor tolerance (which is always in
        ppm) and leave as ``None`` for fragment tolerance (which uses the
        user's ``--tolerance-unit`` preference).

        ``methods_agree`` indicates that multiple independent estimation
        methods (e.g. Gaussian + Crux) converged on the same value.  When
        True, the extreme-divergence guard is bypassed because the consensus
        is strong evidence.

        Cross-unit comparison is always done by converting both values to
        ppm (using ``self.reference_mz``).
        """
        sdrf_val = sdrf_value.get("value") if sdrf_value else None
        sdrf_unit = (sdrf_value.get("unit", "") or "").lower() if sdrf_value else ""

        if detected_value is None:
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=sdrf_value,
                detected_value=None,
                status=ComparisonStatus.MISSING_DETECTED,
            )

        # Determine output unit:
        #   preferred_unit (explicit per-call) > self.tolerance_unit (user flag) > detected_unit
        output_unit = preferred_unit or self.tolerance_unit or detected_unit
        ref_mz = self.reference_mz

        # Convert detected value to output unit
        try:
            output_value = convert_tolerance(detected_value, detected_unit, output_unit, ref_mz)
        except ValueError:
            output_value = detected_value
            output_unit = detected_unit

        # Round for display
        if output_unit.lower() == "ppm":
            output_value = round(output_value, 2)
        else:
            output_value = round(output_value, 4)

        detected_dict = {"value": output_value, "unit": output_unit}

        # ── HEURISTIC-ONLY: never change the SDRF, just recommend review ──
        if tolerance_source == "heuristic":
            heuristic_display = f"{output_value} {output_unit}"
            if sdrf_val is not None:
                recommendation = (
                    f"We recommend checking {parameter_name} because "
                    f"{instrument_name or 'this instrument'} typically uses "
                    f"~{heuristic_display}.  The current SDRF value "
                    f"({sdrf_val} {sdrf_unit}) was kept unchanged because no "
                    f"empirical estimation (Crux/Gaussian) was available."
                )
            else:
                recommendation = (
                    f"We recommend adding {parameter_name} (~{heuristic_display}) "
                    f"based on {instrument_name or 'this instrument'} defaults.  "
                    f"No empirical estimation was available to confirm this value."
                )
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=sdrf_value,
                detected_value=detected_dict,
                status=ComparisonStatus.MISSING_DETECTED,  # treated as "not detected"
                confidence=0.3,
                recommendation=recommendation,
                details={"heuristic_only": True, "tolerance_source": "heuristic"},
            )

        # ── EMPIRICAL: compare properly (convert both to ppm) ──
        if sdrf_val is None:
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=sdrf_value,
                detected_value=detected_dict,
                status=ComparisonStatus.MISSING_SDRF,
                recommendation=f"Add {parameter_name}: {output_value} {output_unit}",
                details={"tolerance_source": "empirical"},
            )

        # Convert both to ppm for comparison
        try:
            sdrf_ppm = convert_tolerance(sdrf_val, sdrf_unit, "ppm", ref_mz) if sdrf_unit else sdrf_val
            detected_ppm = convert_tolerance(detected_value, detected_unit, "ppm", ref_mz)
        except ValueError:
            # Can't convert – fall back to raw numeric comparison
            sdrf_ppm = sdrf_val
            detected_ppm = detected_value

        relative_diff = abs(sdrf_ppm - detected_ppm) / max(sdrf_ppm, detected_ppm, 1e-9)

        # Divergence factor: how many times larger is the bigger value?
        divergence_factor = max(sdrf_ppm, detected_ppm) / max(min(sdrf_ppm, detected_ppm), 1e-9)

        if relative_diff <= TOLERANCE_MATCH_THRESHOLD:
            status = ComparisonStatus.MATCH
            recommendation = None
        elif divergence_factor > TOLERANCE_DIVERGENCE_HARD:
            # Hard limit: >10x divergence is always flagged, even when methods
            # agree.  Such a large gap typically means the SDRF captures a
            # search-engine parameter while the empirical value reflects
            # instrument accuracy — fundamentally different quantities.
            status = ComparisonStatus.MISSING_DETECTED
            recommendation = (
                f"Review {parameter_name}: SDRF has {sdrf_val} {sdrf_unit}, "
                f"empirical estimate is {output_value} {output_unit} "
                f"({divergence_factor:.0f}x divergence). This large difference "
                f"suggests the SDRF value may reflect search-engine settings "
                f"rather than instrument accuracy."
            )
            logger.info(
                f"{parameter_name}: {divergence_factor:.0f}x divergence "
                f"(>{TOLERANCE_DIVERGENCE_HARD:.0f}x hard limit), "
                f"flagging for review"
            )
        elif divergence_factor > TOLERANCE_DIVERGENCE_SOFT and not methods_agree:
            # Soft limit: >5x divergence with weak evidence (single method or
            # methods disagree).  Flag for review.
            status = ComparisonStatus.MISSING_DETECTED
            recommendation = (
                f"Review {parameter_name}: SDRF has {sdrf_val} {sdrf_unit}, "
                f"empirical estimate is {output_value} {output_unit} "
                f"({divergence_factor:.0f}x divergence). Only one estimation "
                f"method available or methods disagree; keeping SDRF unchanged."
            )
            logger.info(
                f"{parameter_name}: {divergence_factor:.0f}x divergence "
                f"with weak consensus, flagging for review"
            )
        elif detected_ppm < sdrf_ppm:
            # Detected is more stringent
            status = ComparisonStatus.IMPROVED
            recommendation = f"Update {parameter_name} to: {output_value} {output_unit}"
        else:
            status = ComparisonStatus.MISMATCH
            recommendation = f"Update {parameter_name} to: {output_value} {output_unit}"

        return ParameterComparison(
            parameter_name=parameter_name,
            sdrf_value=sdrf_value,
            detected_value=detected_dict,
            status=status,
            recommendation=recommendation,
            confidence=0.85,
            details={
                "relative_difference": relative_diff,
                "divergence_factor": divergence_factor,
                "tolerance_source": "empirical",
                "comparison_ppm": {"sdrf": round(sdrf_ppm, 2), "detected": round(detected_ppm, 2)},
            },
        )

    def compare_all(
        self, sdrf_params: Dict[str, Any], detected_params: Dict[str, Any]
    ) -> List[ParameterComparison]:
        """
        Compare all parameters and generate comparison list.

        Args:
            sdrf_params: Parameters extracted from SDRF
            detected_params: Parameters detected from MS analysis

        Returns:
            List of parameter comparisons
        """
        self.comparisons = []

        # Compare instrument
        comp = self.compare_instrument(
            sdrf_params.get("instrument"), detected_params.get("instrument_model")
        )
        self.comparisons.append(comp)

        # Compare MS2 fragmentation/dissociation method
        # Use ms2_fragmentation_type if available, otherwise fall back to fragmentation_type
        ms2_frag_type = detected_params.get("ms2_fragmentation_type") or detected_params.get("fragmentation_type", "unknown")
        comp = self.compare_fragmentation(
            sdrf_params.get("dissociation") or sdrf_params.get("ms2_dissociation"),
            ms2_frag_type,
            param_name="dissociation_method",
        )
        self.comparisons.append(comp)

        # Compare MS3 fragmentation/dissociation method if MS3 data exists
        ms3_frag_type = detected_params.get("ms3_fragmentation_type")
        if ms3_frag_type and ms3_frag_type != "unknown":
            comp = self.compare_fragmentation(
                sdrf_params.get("ms3_dissociation"),
                ms3_frag_type,
                param_name="ms3_dissociation_method",
            )
            self.comparisons.append(comp)

        # Instrument name for heuristic recommendation messages
        inst_model = detected_params.get("instrument_model", {})
        instrument_name = inst_model.get("name", "") if isinstance(inst_model, dict) else str(inst_model)

        # Check whether multiple tolerance estimation methods agreed.
        # The ConsensusResult.recommendation_basis contains "all agree" or
        # "agreeing methods" when Gaussian + Crux converge.
        tol_estimation = detected_params.get("tolerance_estimation")
        tol_basis = getattr(tol_estimation, "recommendation_basis", "") or ""
        prec_methods_agree = "all agree" in tol_basis or "agreeing methods" in tol_basis
        # Fragment basis is after "; fragment: " if present
        frag_basis_part = tol_basis.split("; fragment: ")[-1] if "; fragment: " in tol_basis else tol_basis
        frag_methods_agree = "all agree" in frag_basis_part or "agreeing methods" in frag_basis_part

        # Compare precursor tolerance
        # Precursor is ALWAYS in ppm – the user's --tolerance-unit does NOT apply.
        prec_source = detected_params.get("precursor_tolerance_source", "heuristic")
        comp = self.compare_tolerance(
            "precursor_mass_tolerance",
            sdrf_params.get("precursor_tolerance"),
            detected_params.get("precursor_tolerance_ppm"),
            "ppm",
            tolerance_source=prec_source,
            instrument_name=instrument_name,
            preferred_unit="ppm",  # always ppm for precursor
            methods_agree=prec_methods_agree,
        )
        self.comparisons.append(comp)

        # Compare fragment tolerance
        # Fragment uses the user's --tolerance-unit preference (preferred_unit=None
        # → falls through to self.tolerance_unit → detected_unit).
        frag_source = detected_params.get("fragment_tolerance_source", "heuristic")
        if detected_params.get("fragment_tolerance_ppm") is not None:
            comp = self.compare_tolerance(
                "fragment_mass_tolerance",
                sdrf_params.get("fragment_tolerance"),
                detected_params.get("fragment_tolerance_ppm"),
                "ppm",
                tolerance_source=frag_source,
                instrument_name=instrument_name,
                methods_agree=frag_methods_agree,
                # preferred_unit=None → uses self.tolerance_unit if set
            )
        elif detected_params.get("fragment_tolerance_da") is not None:
            comp = self.compare_tolerance(
                "fragment_mass_tolerance",
                sdrf_params.get("fragment_tolerance"),
                detected_params.get("fragment_tolerance_da"),
                "Da",
                tolerance_source=frag_source,
                instrument_name=instrument_name,
                methods_agree=frag_methods_agree,
            )
        else:
            comp = ParameterComparison(
                parameter_name="fragment_mass_tolerance",
                sdrf_value=sdrf_params.get("fragment_tolerance"),
                detected_value=None,
                status=ComparisonStatus.MISSING_DETECTED,
            )
        self.comparisons.append(comp)

        # Compare MS min charge
        comp = self.compare_charge(
            "ms_min_charge",
            sdrf_params.get("ms_min_charge"),
            detected_params.get("ms_min_charge"),
        )
        self.comparisons.append(comp)

        # Compare MS max charge
        comp = self.compare_charge(
            "ms_max_charge",
            sdrf_params.get("ms_max_charge"),
            detected_params.get("ms_max_charge"),
        )
        self.comparisons.append(comp)

        # Compare collision energy
        comp = self.compare_collision_energy(
            sdrf_params.get("collision_energy"),
            detected_params.get("collision_energy"),
        )
        self.comparisons.append(comp)

        # Compare MS2 mass analyzer
        comp = self.compare_ontology_value(
            "ms2_mass_analyzer",
            sdrf_params.get("ms2_mass_analyzer"),
            detected_params.get("mass_analyzer_type"),
        )
        self.comparisons.append(comp)

        # Compare scan window lower/upper limits (from ms1_scan_range tuple)
        ms1_scan_range = detected_params.get("ms1_scan_range")
        if ms1_scan_range and isinstance(ms1_scan_range, (list, tuple)) and len(ms1_scan_range) >= 2:
            comp = self.compare_numeric_mz(
                "scan_window_lower",
                sdrf_params.get("scan_window_lower"),
                float(ms1_scan_range[0]),
            )
            self.comparisons.append(comp)

            comp = self.compare_numeric_mz(
                "scan_window_upper",
                sdrf_params.get("scan_window_upper"),
                float(ms1_scan_range[1]),
            )
            self.comparisons.append(comp)

        # Compare isolation window width
        isolation_window = detected_params.get("isolation_window_da")
        if isolation_window is not None:
            comp = self.compare_numeric_mz(
                "isolation_window_width",
                sdrf_params.get("isolation_window_width"),
                float(isolation_window),
            )
            self.comparisons.append(comp)

        return self.comparisons

    def compare_collision_energy(
        self,
        sdrf_value: Optional[str],
        detected_value: Optional[str],
    ) -> ParameterComparison:
        """
        Compare collision energy values.

        SDRF spec (dda-acquisition, dia-acquisition):
        Pattern: ^\\d+(\\.\\d+)?%? (NCE|eV)(;\\d+(\\.\\d+)?%? (NCE|eV))*$|^not available$
        Valid: "30 NCE", "30% NCE", "27 eV", "25 NCE;27 NCE;30 NCE", "not available"
        """
        if not detected_value or not detected_value.strip():
            return ParameterComparison(
                parameter_name="collision_energy",
                sdrf_value=sdrf_value,
                detected_value=None,
                status=ComparisonStatus.MISSING_DETECTED,
            )

        sdrf_val = (sdrf_value or "").strip()
        # Treat "not available" / "not applicable" as empty (spec allows these)
        if sdrf_val.lower() in ("not available", "not applicable"):
            sdrf_val = ""

        if not sdrf_val:
            return ParameterComparison(
                parameter_name="collision_energy",
                sdrf_value=sdrf_value,
                detected_value=detected_value,
                status=ComparisonStatus.MISSING_SDRF,
                recommendation=f"Add collision energy: {detected_value}",
            )

        # Exact match
        if sdrf_val.lower() == detected_value.strip().lower():
            return ParameterComparison(
                parameter_name="collision_energy",
                sdrf_value=sdrf_value,
                detected_value=detected_value,
                status=ComparisonStatus.MATCH,
            )

        # Same numeric value, different unit (e.g. 32 NCE vs 32 eV) - treat as match
        sdrf_match = re.match(r"^([\d.]+)\s*(%?\s*(?:NCE|eV))", sdrf_val, re.I)
        det_match = re.match(r"^([\d.]+)\s*(%?\s*(?:NCE|eV))", detected_value, re.I)
        if sdrf_match and det_match:
            sdrf_num = float(sdrf_match.group(1))
            det_num = float(det_match.group(1))
            if abs(sdrf_num - det_num) < 0.01:
                return ParameterComparison(
                    parameter_name="collision_energy",
                    sdrf_value=sdrf_value,
                    detected_value=detected_value,
                    status=ComparisonStatus.MATCH,
                    details={"unit_mismatch": True},
                )

        return ParameterComparison(
            parameter_name="collision_energy",
            sdrf_value=sdrf_value,
            detected_value=detected_value,
            status=ComparisonStatus.MISMATCH,
            recommendation=f"Update collision energy to: {detected_value}",
        )

    def compare_ontology_value(
        self,
        parameter_name: str,
        sdrf_value: Optional[Dict],
        detected_value: Optional[Dict],
    ) -> ParameterComparison:
        """
        Compare a generic ontology-valued parameter (e.g. mass analyzer).

        Both ``sdrf_value`` and ``detected_value`` are expected to be dicts
        with ``name`` and optionally ``accession`` keys.
        """
        detected_name = ""
        if isinstance(detected_value, dict):
            detected_name = (detected_value.get("name", "") or "").lower()
        elif isinstance(detected_value, str):
            detected_name = detected_value.lower()
            detected_value = {"name": detected_value}

        if not detected_name or detected_name == "unknown":
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=sdrf_value,
                detected_value=detected_value,
                status=ComparisonStatus.MISSING_DETECTED,
            )

        sdrf_name = ""
        if isinstance(sdrf_value, dict):
            sdrf_name = (sdrf_value.get("name", "") or "").lower()
        elif isinstance(sdrf_value, str):
            sdrf_name = sdrf_value.lower()

        if not sdrf_name:
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=sdrf_value,
                detected_value=detected_value,
                status=ComparisonStatus.MISSING_SDRF,
                recommendation=(
                    f"Add {parameter_name}: "
                    f"NT={detected_value.get('name', '')};AC={detected_value.get('accession', '')}"
                ),
            )

        # Fuzzy matching (contains or equals)
        if sdrf_name == detected_name or sdrf_name in detected_name or detected_name in sdrf_name:
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=sdrf_value,
                detected_value=detected_value,
                status=ComparisonStatus.MATCH,
            )

        return ParameterComparison(
            parameter_name=parameter_name,
            sdrf_value=sdrf_value,
            detected_value=detected_value,
            status=ComparisonStatus.MISMATCH,
            recommendation=(
                f"Update {parameter_name} to: "
                f"NT={detected_value.get('name', '')};AC={detected_value.get('accession', '')}"
            ),
        )

    def compare_charge(
        self,
        parameter_name: str,
        sdrf_value: Optional[int],
        detected_value: Optional[int],
    ) -> ParameterComparison:
        """Compare charge state values (min or max)."""
        if detected_value is None:
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=sdrf_value,
                detected_value=None,
                status=ComparisonStatus.MISSING_DETECTED,
            )

        if sdrf_value is None:
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=None,
                detected_value=detected_value,
                status=ComparisonStatus.MISSING_SDRF,
                recommendation=f"Add {parameter_name}: {detected_value}",
            )

        if sdrf_value == detected_value:
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=sdrf_value,
                detected_value=detected_value,
                status=ComparisonStatus.MATCH,
            )
        else:
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=sdrf_value,
                detected_value=detected_value,
                status=ComparisonStatus.MISMATCH,
                recommendation=f"Update {parameter_name} to: {detected_value}",
            )

    def compare_numeric_mz(
        self,
        parameter_name: str,
        sdrf_value: Optional[float],
        detected_value: Optional[float],
    ) -> ParameterComparison:
        """Compare simple numeric m/z values (scan windows, isolation width)."""
        if detected_value is None:
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=sdrf_value,
                detected_value=None,
                status=ComparisonStatus.MISSING_DETECTED,
            )

        if sdrf_value is None:
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=None,
                detected_value=detected_value,
                status=ComparisonStatus.MISSING_SDRF,
                recommendation=f"Add {parameter_name}: {detected_value:g}",
            )

        # Values within 1% are considered a match
        relative_diff = abs(sdrf_value - detected_value) / max(abs(sdrf_value), abs(detected_value), 1e-9)
        if relative_diff <= 0.01:
            return ParameterComparison(
                parameter_name=parameter_name,
                sdrf_value=sdrf_value,
                detected_value=detected_value,
                status=ComparisonStatus.MATCH,
            )

        return ParameterComparison(
            parameter_name=parameter_name,
            sdrf_value=sdrf_value,
            detected_value=detected_value,
            status=ComparisonStatus.MISMATCH,
            recommendation=f"Update {parameter_name} to: {detected_value:g}",
        )

    def get_refinements(self) -> List[ParameterComparison]:
        """Get list of parameters that should be refined."""
        return [
            c
            for c in self.comparisons
            if c.status
            in [
                ComparisonStatus.MISMATCH,
                ComparisonStatus.MISSING_SDRF,
                ComparisonStatus.IMPROVED,
            ]
        ]
