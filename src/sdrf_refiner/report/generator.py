"""
ReportGenerator - Generate plain text refinement reports.
"""

import logging
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from sdrf_refiner.sdrf.comparer import ComparisonStatus, ParameterComparison

logger = logging.getLogger(__name__)


class ReportGenerator:
    """
    Generates detailed plain text refinement reports showing:
    - Detected vs original parameters
    - Confidence levels for each recommendation
    - Supporting evidence from analysis
    """

    SEPARATOR = "=" * 60

    def __init__(self):
        self.report_lines: List[str] = []

    def generate(
        self,
        pride_accession: str,
        sdrf_params: Dict[str, Any],
        detected_params: Dict[str, Any],
        comparisons: List[ParameterComparison],
        analysis_results: List,
        output_path: Path,
        changes: Optional[List[Dict[str, Any]]] = None,
        ptm_results: Optional[List] = None,
    ) -> str:
        """
        Generate comprehensive refinement report.

        Args:
            pride_accession: PRIDE accession number
            sdrf_params: Original parameters from SDRF
            detected_params: Parameters detected from analysis
            comparisons: List of parameter comparisons
            analysis_results: List of per-file analysis results
            output_path: Path to write report
            changes: List of changes applied to SDRF
            ptm_results: Per-file PTM detection results (optional)

        Returns:
            Report content as string
        """
        self.report_lines = []
        self._changes = changes or []
        self._detected_params = detected_params
        self._sdrf_params = sdrf_params
        self._ptm_results = ptm_results or []

        # Header
        self._add_header(pride_accession, len(analysis_results))

        # Summary
        self._add_summary(comparisons)

        # Original SDRF parameters
        self._add_original_params(sdrf_params)

        # Detected parameters
        self._add_detected_params(detected_params)

        # Tolerance estimation details
        self._add_tolerance_estimation(detected_params)

        # PTM detection results
        self._add_ptm_detection_section(detected_params)

        # Changes made (new section with column: old -> new format)
        self._add_changes_section()

        # Recommendations
        self._add_recommendations(comparisons)

        # Detailed file results
        self._add_file_results(analysis_results)

        # Footer
        self._add_footer()

        # Join and write
        report_content = "\n".join(self.report_lines)

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        with open(output_path, "w", encoding="utf-8") as f:
            f.write(report_content)

        logger.info(f"Report written to: {output_path}")

        # Also log summary to console
        self._log_console_summary(comparisons)

        return report_content

    def _add_line(self, line: str = "") -> None:
        """Add a line to the report."""
        self.report_lines.append(line)

    def _add_header(self, pride_accession: str, num_files: int) -> None:
        """Add report header."""
        self._add_line(self.SEPARATOR)
        self._add_line("SDRF REFINEMENT REPORT")
        self._add_line(self.SEPARATOR)
        self._add_line(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self._add_line(f"Accession: {pride_accession}")
        self._add_line(f"Files analyzed: {num_files}")
        self._add_line(self.SEPARATOR)
        self._add_line()

    def _add_summary(self, comparisons: List[ParameterComparison]) -> None:
        """Add summary section."""
        status_counts = {}
        for comp in comparisons:
            status = comp.status.value
            status_counts[status] = status_counts.get(status, 0) + 1

        self._add_line("SUMMARY")
        self._add_line("-" * 40)
        self._add_line(f"Parameters compared: {len(comparisons)}")
        self._add_line(f"  Matching: {status_counts.get('match', 0)}")
        self._add_line(f"  Mismatches: {status_counts.get('mismatch', 0)}")
        self._add_line(f"  Missing in SDRF: {status_counts.get('missing_in_sdrf', 0)}")
        self._add_line(f"  Can improve: {status_counts.get('can_improve', 0)}")
        self._add_line(f"  Could not detect: {status_counts.get('missing_detected', 0)}")

        needs_refinement = (
            status_counts.get("mismatch", 0)
            + status_counts.get("missing_in_sdrf", 0)
            + status_counts.get("can_improve", 0)
        )
        if needs_refinement > 0:
            self._add_line()
            self._add_line(f"*** {needs_refinement} parameter(s) recommended for refinement ***")
        else:
            self._add_line()
            self._add_line("No refinements needed - SDRF matches detected parameters.")

        self._add_line()

    def _add_original_params(self, sdrf_params: Dict[str, Any]) -> None:
        """Add original SDRF parameters section."""
        self._add_line("ORIGINAL SDRF PARAMETERS")
        self._add_line("-" * 40)

        if not sdrf_params:
            self._add_line("  (no parameters found in SDRF)")
            self._add_line()
            return

        for key, value in sdrf_params.items():
            if value:
                formatted = self._format_param_value(value)
                self._add_line(f"  {key}: {formatted}")
            else:
                self._add_line(f"  {key}: (not specified)")

        self._add_line()

    def _add_detected_params(self, detected_params: Dict[str, Any]) -> None:
        """Add detected parameters section."""
        self._add_line("DETECTED PARAMETERS")
        self._add_line("-" * 40)

        if not detected_params:
            self._add_line("  (no parameters detected)")
            self._add_line()
            return

        # Instrument
        inst = detected_params.get("instrument_model")
        if isinstance(inst, dict):
            name = inst.get("name", "unknown")
            accession = inst.get("accession", "")
            self._add_line(f"  Instrument: {name} ({accession})")
        elif inst:
            self._add_line(f"  Instrument: {inst}")

        # Fragmentation
        ms2_frag = detected_params.get("ms2_fragmentation_type") or detected_params.get("fragmentation_type")
        ms3_frag = detected_params.get("ms3_fragmentation_type")
        has_ms3 = ms3_frag and ms3_frag != "unknown"

        if ms2_frag:
            # Only add "MS2" prefix if MS3 is also present
            if has_ms3:
                self._add_line(f"  MS2 Fragmentation: {ms2_frag}")
            else:
                self._add_line(f"  Fragmentation: {ms2_frag}")

        if has_ms3:
            self._add_line(f"  MS3 Fragmentation: {ms3_frag}")

        # Tolerances (with source annotation)
        prec_tol = detected_params.get("precursor_tolerance_ppm")
        prec_src = detected_params.get("precursor_tolerance_source", "")
        if prec_tol is not None:
            src_label = f" [{prec_src}]" if prec_src else ""
            self._add_line(f"  Precursor tolerance: {prec_tol:.1f} ppm{src_label}")

        frag_tol_ppm = detected_params.get("fragment_tolerance_ppm")
        frag_tol_da = detected_params.get("fragment_tolerance_da")
        frag_src = detected_params.get("fragment_tolerance_source", "")
        src_label = f" [{frag_src}]" if frag_src else ""
        if frag_tol_ppm is not None:
            self._add_line(f"  Fragment tolerance: {frag_tol_ppm:.1f} ppm{src_label}")
        elif frag_tol_da is not None:
            self._add_line(f"  Fragment tolerance: {frag_tol_da:.4f} Da{src_label}")

        # Acquisition type
        acq = detected_params.get("acquisition_type")
        if acq:
            self._add_line(f"  Acquisition type: {acq}")

        # Charge state range
        min_charge = detected_params.get("ms_min_charge")
        max_charge = detected_params.get("ms_max_charge")
        if min_charge is not None and max_charge is not None:
            self._add_line(f"  Charge state range: {min_charge}-{max_charge}")
        elif min_charge is not None:
            self._add_line(f"  Min charge: {min_charge}")
        elif max_charge is not None:
            self._add_line(f"  Max charge: {max_charge}")

        # Collision energy
        ce = detected_params.get("collision_energy")
        if ce:
            self._add_line(f"  Collision energy: {ce}")

        # pyopenms-derived properties
        ma = detected_params.get("mass_analyzer_type")
        if ma:
            if isinstance(ma, dict):
                name = ma.get("name", "")
                acc = ma.get("accession", "")
                self._add_line(f"  Mass analyzer: {name} ({acc})" if acc else f"  Mass analyzer: {name}")
            elif ma != "unknown":
                self._add_line(f"  Mass analyzer: {ma}")

        ionization = detected_params.get("ionization_method")
        if ionization:
            if isinstance(ionization, dict):
                name = ionization.get("name", "")
                acc = ionization.get("accession", "")
                self._add_line(f"  Ionization: {name} ({acc})" if acc else f"  Ionization: {name}")
            elif ionization != "unknown":
                self._add_line(f"  Ionization: {ionization}")

        polarity = detected_params.get("polarity")
        if polarity and polarity != "unknown":
            self._add_line(f"  Polarity: {polarity}")

        scan_range = detected_params.get("ms1_scan_range")
        if scan_range:
            self._add_line(f"  MS1 scan range: {scan_range[0]:.0f}-{scan_range[1]:.0f} m/z")

        isol_win = detected_params.get("isolation_window_da")
        if isol_win is not None:
            self._add_line(f"  Isolation window: {isol_win:.2f} Da")

        serial = detected_params.get("instrument_serial_number")
        if serial:
            self._add_line(f"  Instrument serial: {serial}")

        self._add_line()

    def _add_tolerance_estimation(self, detected_params: Dict[str, Any]) -> None:
        """Add tolerance estimation details section."""
        tolerance_est = detected_params.get("tolerance_estimation")
        prec_source = detected_params.get("precursor_tolerance_source", "unknown")
        frag_source = detected_params.get("fragment_tolerance_source", "unknown")

        self._add_line("TOLERANCE ESTIMATION")
        self._add_line("-" * 40)

        if tolerance_est and hasattr(tolerance_est, "format_summary"):
            ref_mz = detected_params.get("reference_mz", 500.0)
            for line in tolerance_est.format_summary(ref_mz).splitlines():
                self._add_line(f"  {line}")
        else:
            self._add_line("  No tolerance estimation details available.")

        self._add_line(f"  Precursor tolerance source: {prec_source}")
        self._add_line(f"  Fragment tolerance source: {frag_source}")

        if prec_source == "heuristic" or frag_source == "heuristic":
            self._add_line()
            self._add_line(
                "  NOTE: Heuristic tolerances are based on instrument category "
                "defaults and have NOT been empirically verified.  SDRF values "
                "were kept unchanged for heuristic-only parameters."
            )

        self._add_line()

    def _build_sdrf_unimod_set(self) -> set:
        """Build a set of uppercased UNIMOD accessions from SDRF modifications."""
        sdrf_mods = self._sdrf_params.get("modifications", [])
        return {
            m["accession"].upper()
            for m in sdrf_mods
            if m.get("accession")
        }

    def _add_ptm_detection_section(self, detected_params: Dict[str, Any]) -> None:
        """Add PTM detection results section with SDRF comparison."""
        ptm_data = detected_params.get("ptm_detection")
        if not ptm_data:
            return

        hits = ptm_data.get("hits", [])
        if not hits:
            return

        n_files = ptm_data.get("n_files", 0)
        sdrf_unimod_set = self._build_sdrf_unimod_set()

        self._add_line("PTM DETECTION")
        self._add_line("-" * 40)
        self._add_line(f"  Files analyzed for PTMs: {n_files}")
        self._add_line()

        # Collect all detected UNIMOD accessions for the "not detected" section
        detected_unimod_set: set = set()

        tier_names = {1: "Reporter Ions", 2: "Diagnostic Ions", 3: "Mass Shifts"}
        for tier in [1, 2, 3]:
            tier_hits = [h for h in hits if h["tier"] == tier]
            if not tier_hits:
                continue

            self._add_line(f"  {tier_names[tier]}:")
            for hit in tier_hits:
                unimod = hit.get("unimod_accession", "")
                if unimod:
                    detected_unimod_set.add(unimod.upper())
                self._add_line(
                    f"    {self._format_ptm_hit(hit, tier, n_files, sdrf_unimod_set)}"
                )
            self._add_line()

        # SDRF-only modifications (not detected from spectra)
        sdrf_mods = self._sdrf_params.get("modifications", [])
        sdrf_only = [
            m for m in sdrf_mods
            if m.get("accession") and m["accession"].upper() not in detected_unimod_set
        ]
        if sdrf_only:
            self._add_line("  Not detected (in SDRF only):")
            for m in sdrf_only:
                name = m.get("name", "unknown")
                acc = m.get("accession", "")
                mt = m.get("modification_type", "")
                mt_str = f" [{mt}]" if mt else ""
                self._add_line(f"    {name} ({acc}){mt_str}")
            self._add_line()

        # Per-file PTM summary
        if self._ptm_results:
            self._add_line("  Per-file PTM Summary:")
            for pr in self._ptm_results:
                fname = getattr(pr.mzml_file, "name", str(pr.mzml_file))
                if pr.hits:
                    parts = [
                        f"{h.name} ({h.evidence_count}/{pr.n_spectra_scanned})"
                        for h in pr.hits
                    ]
                    self._add_line(f"    {fname}: {', '.join(parts)}")
                else:
                    self._add_line(f"    {fname}: no PTMs detected")
            self._add_line()

    def _format_ptm_hit(
        self,
        hit: Dict[str, Any],
        tier: int,
        n_files: int,
        sdrf_unimod_set: set,
    ) -> str:
        """Format a single PTM hit for the report."""
        name = hit["name"]
        unimod = hit["unimod_accession"]
        evidence = hit["total_evidence"]
        scanned = hit["total_scanned"]
        conf_label = self._confidence_label(hit["max_confidence"])
        files_det = hit["files_detected"]
        mod_type = hit.get("modification_type", "variable")

        # Compare against SDRF
        sdrf_status = "match" if unimod.upper() in sdrf_unimod_set else "new"

        if tier == 1:
            pct = (evidence / scanned * 100) if scanned > 0 else 0
            return (
                f"{name} ({unimod}): detected in {pct:.0f}% of spectra "
                f"({files_det}/{n_files} files) [{conf_label}] [{mod_type}] [{sdrf_status}]"
            )

        if tier == 2:
            pct = (evidence / scanned * 100) if scanned > 0 else 0
            return (
                f"{name} ({unimod}): {pct:.1f}% of spectra "
                f"({evidence} observations) [{conf_label}] [{mod_type}] [{sdrf_status}]"
            )

        # Tier 3: mass shifts with enrichment stats
        delta = hit.get("mass_delta")
        delta_str = f" ({delta:+.3f} Da)" if delta is not None else ""

        stats_parts = []
        enrich = hit.get("max_enrichment")
        prob = hit.get("max_probability")
        if enrich is not None:
            stats_parts.append(f"enrichment={enrich:.1f}x")
        if prob is not None:
            stats_parts.append(f"prob={prob:.2f}")
        stats_str = f" ({', '.join(stats_parts)})" if stats_parts else ""

        return (
            f"{name}{delta_str} ({unimod}): {evidence} observations "
            f"({files_det}/{n_files} files) [{conf_label}] [{mod_type}] [{sdrf_status}]{stats_str}"
        )

    def _add_changes_section(self) -> None:
        """Add section showing changes made in column: old -> new format."""
        self._add_line("CHANGES APPLIED")
        self._add_line("-" * 40)

        if not self._changes:
            self._add_line("  No changes applied.")
            self._add_line()
            return

        for change in self._changes:
            col = change.get("column", "")
            old_val = change.get("old_value", "(not specified)")
            new_val = change.get("new_value", "")
            self._add_line(f"  {col}: {old_val} -> {new_val}")

        self._add_line()

    def _add_recommendations(self, comparisons: List[ParameterComparison]) -> None:
        """Add recommendations section."""
        self._add_line("RECOMMENDATIONS")
        self._add_line("-" * 40)

        # Active refinements (will be applied to SDRF)
        refinements = [
            c
            for c in comparisons
            if c.status
            in [
                ComparisonStatus.MISMATCH,
                ComparisonStatus.MISSING_SDRF,
                ComparisonStatus.IMPROVED,
            ]
            and not c.details.get("heuristic_only")
        ]

        # Heuristic-only review items (will NOT be applied, only shown)
        heuristic_reviews = [
            c for c in comparisons if c.details.get("heuristic_only")
        ]

        if not refinements and not heuristic_reviews:
            self._add_line("  No refinements recommended.")
            self._add_line()
            return

        # Sort by priority
        refinements.sort(key=lambda x: self._get_priority(x))

        for ref in refinements:
            confidence = self._confidence_label(ref.confidence)
            col_name = self._param_to_column(ref.parameter_name)

            # Format as: comment[column]: old -> new
            if ref.sdrf_value:
                old_val = self._format_param_value(ref.sdrf_value)
            else:
                old_val = "(not specified)"

            if ref.detected_value:
                new_val = self._format_param_value(ref.detected_value)
            else:
                new_val = "(unknown)"

            self._add_line(f"[{confidence}] {col_name}: {old_val} -> {new_val}")

        # Heuristic-only items get a separate [REVIEW] block
        if heuristic_reviews:
            self._add_line()
            self._add_line("REVIEW (heuristic only - SDRF values kept unchanged):")
            for ref in heuristic_reviews:
                if ref.recommendation:
                    col_name = self._param_to_column(ref.parameter_name)
                    self._add_line(f"  [REVIEW] {col_name}: {ref.recommendation}")

        self._add_line()

    def _param_to_column(self, param_name: str) -> str:
        """Map parameter name to SDRF column name."""
        from sdrf_refiner.config import SDRF_COLUMNS
        mapping = {
            "instrument": SDRF_COLUMNS["instrument"],
            "dissociation_method": SDRF_COLUMNS["dissociation"],
            "ms3_dissociation_method": SDRF_COLUMNS["ms3_dissociation"],
            "precursor_mass_tolerance": SDRF_COLUMNS["precursor_tolerance"],
            "fragment_mass_tolerance": SDRF_COLUMNS["fragment_tolerance"],
            "ms_min_charge": SDRF_COLUMNS["ms_min_charge"],
            "ms_max_charge": SDRF_COLUMNS["ms_max_charge"],
            "collision_energy": SDRF_COLUMNS["collision_energy"],
            "ms2_mass_analyzer": SDRF_COLUMNS["ms2_mass_analyzer"],
        }
        return mapping.get(param_name, param_name)

    def _add_file_results(self, analysis_results: List) -> None:
        """Add detailed per-file results section."""
        if not analysis_results:
            return

        self._add_line("DETAILED FILE RESULTS")
        self._add_line("-" * 40)

        for result in analysis_results:
            filename = result.mzml_file.name if hasattr(result.mzml_file, "name") else str(result.mzml_file)
            self._add_line(f"File: {filename}")

            inst = result.instrument_model
            if inst:
                name = inst.get("name", "unknown") if isinstance(inst, dict) else str(inst)
                self._add_line(f"  Instrument: {name}")

            # Fragmentation
            ms2_frag = getattr(result, 'ms2_fragmentation_type', None) or result.fragmentation_type
            ms3_frag = getattr(result, 'ms3_fragmentation_type', None)
            has_ms3 = ms3_frag and ms3_frag != "unknown"

            if has_ms3:
                self._add_line(f"  MS2 Fragmentation: {ms2_frag}")
                n_ms3 = getattr(result, 'n_ms3_spectra', 0)
                self._add_line(f"  MS3 Fragmentation: {ms3_frag} ({n_ms3} spectra)")
            else:
                self._add_line(f"  Fragmentation: {ms2_frag}")

            self._add_line(f"  High accuracy: {'Yes' if result.high_accuracy_precursors else 'No'}")

            if result.precursor_tolerance_ppm:
                self._add_line(f"  Precursor tol: {result.precursor_tolerance_ppm:.1f} ppm")

            if result.fragment_tolerance_ppm:
                self._add_line(f"  Fragment tol: {result.fragment_tolerance_ppm:.1f} ppm")
            elif result.fragment_tolerance_da:
                self._add_line(f"  Fragment tol: {result.fragment_tolerance_da:.4f} Da")

            # Charge state range
            min_charge = getattr(result, 'ms_min_charge', None)
            max_charge = getattr(result, 'ms_max_charge', None)
            if min_charge is not None and max_charge is not None:
                self._add_line(f"  Charge range: {min_charge}-{max_charge}")

            # Collision energy
            ce = getattr(result, 'collision_energy', None)
            if ce:
                self._add_line(f"  Collision energy: {ce}")

            # pyopenms-derived properties
            ma = getattr(result, 'mass_analyzer_type', 'unknown')
            if ma != "unknown":
                self._add_line(f"  Mass analyzer: {ma}")

            ionization = getattr(result, 'ionization_method', 'unknown')
            if ionization != "unknown":
                self._add_line(f"  Ionization: {ionization}")

            polarity = getattr(result, 'polarity', 'unknown')
            if polarity != "unknown":
                self._add_line(f"  Polarity: {polarity}")

            scan_range = getattr(result, 'ms1_scan_range', None)
            if scan_range:
                self._add_line(f"  MS1 scan range: {scan_range[0]:.0f}-{scan_range[1]:.0f} m/z")

            isol_win = getattr(result, 'isolation_window_da', None)
            if isol_win is not None:
                self._add_line(f"  Isolation window: {isol_win:.2f} Da")

            serial = getattr(result, 'instrument_serial_number', None)
            if serial:
                self._add_line(f"  Instrument serial: {serial}")

            n_ms3 = getattr(result, 'n_ms3_spectra', 0)
            if n_ms3 > 0:
                self._add_line(f"  Spectra: {result.n_spectra} (MS1: {result.n_ms1_spectra}, MS2: {result.n_ms2_spectra}, MS3: {n_ms3})")
            else:
                self._add_line(f"  Spectra: {result.n_spectra} (MS1: {result.n_ms1_spectra}, MS2: {result.n_ms2_spectra})")
            self._add_line()

    def _add_footer(self) -> None:
        """Add report footer."""
        self._add_line(self.SEPARATOR)
        self._add_line("Generated by sdrf-refiner v1.0.0")
        self._add_line(self.SEPARATOR)

    def _format_param_value(self, value: Any) -> str:
        """Format a parameter value for display."""
        if value is None:
            return "(not specified)"

        if isinstance(value, dict):
            # Ontology value
            name = value.get("name") or value.get("raw", "")
            accession = value.get("accession", "")
            unit = value.get("unit", "")
            val = value.get("value")

            if val is not None and unit:
                return f"{val} {unit}"
            elif name and accession:
                return f"{name} ({accession})"
            elif name:
                return name
            else:
                return str(value)

        return str(value)

    def _confidence_label(self, score: float) -> str:
        """Convert confidence score to human-readable label."""
        if score >= 0.9:
            return "HIGH"
        elif score >= 0.7:
            return "MEDIUM"
        else:
            return "LOW"

    def _get_priority(self, comp: ParameterComparison) -> int:
        """Assign priority to refinement (lower = higher priority)."""
        if comp.status == ComparisonStatus.MISMATCH:
            return 1
        elif comp.status == ComparisonStatus.MISSING_SDRF:
            return 2
        else:
            return 3

    def _log_console_summary(self, comparisons: List[ParameterComparison]) -> None:
        """Log a brief summary to console."""
        refinements = [
            c
            for c in comparisons
            if c.status
            in [
                ComparisonStatus.MISMATCH,
                ComparisonStatus.MISSING_SDRF,
                ComparisonStatus.IMPROVED,
            ]
            and not c.details.get("heuristic_only")
        ]

        heuristic_reviews = [
            c for c in comparisons if c.details.get("heuristic_only")
        ]

        logger.info(self.SEPARATOR)
        logger.info("SDRF REFINEMENT SUMMARY")
        logger.info(self.SEPARATOR)

        status_counts = {}
        for comp in comparisons:
            if comp.details.get("heuristic_only"):
                continue  # Don't count heuristic-only in normal stats
            status = comp.status.value
            status_counts[status] = status_counts.get(status, 0) + 1

        logger.info(f"Parameters compared: {len(comparisons)}")
        logger.info(f"  - Matching: {status_counts.get('match', 0)}")
        logger.info(f"  - Mismatches: {status_counts.get('mismatch', 0)}")
        logger.info(f"  - Missing in SDRF: {status_counts.get('missing_in_sdrf', 0)}")
        logger.info(f"  - Can improve: {status_counts.get('can_improve', 0)}")
        if heuristic_reviews:
            logger.info(f"  - Heuristic only (review): {len(heuristic_reviews)}")

        # Log changes in column: old -> new format
        if self._changes:
            logger.info("")
            logger.info("CHANGES APPLIED:")
            for change in self._changes:
                col = change.get("column", "")
                old_val = change.get("old_value", "(not specified)")
                new_val = change.get("new_value", "")
                logger.info(f"  {col}: {old_val} -> {new_val}")
        elif refinements:
            logger.info("")
            logger.info("RECOMMENDATIONS:")
            for ref in refinements:
                confidence = self._confidence_label(ref.confidence)
                col_name = self._param_to_column(ref.parameter_name)
                old_val = self._format_param_value(ref.sdrf_value) if ref.sdrf_value else "(not specified)"
                new_val = self._format_param_value(ref.detected_value) if ref.detected_value else "(unknown)"
                logger.info(f"  [{confidence}] {col_name}: {old_val} -> {new_val}")
        else:
            logger.info("")
            logger.info("No refinements recommended - SDRF matches detected parameters.")

        # Log heuristic reviews
        if heuristic_reviews:
            logger.info("")
            logger.info("REVIEW (heuristic only - no changes applied):")
            for ref in heuristic_reviews:
                col_name = self._param_to_column(ref.parameter_name)
                if ref.recommendation:
                    logger.info(f"  [REVIEW] {ref.recommendation}")

        logger.info(self.SEPARATOR)
