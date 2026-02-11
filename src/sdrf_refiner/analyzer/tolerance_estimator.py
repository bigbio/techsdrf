"""
Tolerance Estimator - Empirical mass tolerance estimation using multiple methods.

Integrates two approaches:
1. Gaussian fitting (adapted from RunAssessor) - fits mass error distribution
2. Crux param-medic - external tool for tolerance prediction

Provides consensus recommendation when methods disagree.
"""

import logging
import os
import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# =============================================================================
# Unit Conversion Utilities
# =============================================================================

# Canonical unit names and aliases (Da-like: da, Da, th, Th, m/z)
_DA_ALIASES = frozenset({"da", "dalton", "daltons", "th", "thompson", "m/z"})
_PPM_ALIASES = frozenset({"ppm", "parts per million"})


def _normalize_unit(unit: str) -> str:
    """Normalize unit string to canonical form ('ppm' or 'Da')."""
    u = str(unit).strip().lower()
    if u in _PPM_ALIASES or u == "ppm":
        return "ppm"
    if u in _DA_ALIASES or u == "da":
        return "Da"
    raise ValueError(f"Unknown mass tolerance unit: {unit!r}. Use 'ppm' or 'Da' (aliases: th, m/z).")


def convert_tolerance(
    value: float,
    from_unit: str,
    to_unit: str,
    reference_mz: float = 500.0,
) -> float:
    """
    Convert a mass tolerance value between ppm and Da.

    Use this when an algorithm provides a value in one unit and you need it
    in another. Conversion requires a reference m/z (ppm is relative).

    Args:
        value: Tolerance value in from_unit
        from_unit: Source unit ('ppm', 'Da', 'Th', 'm/z', etc.)
        to_unit: Target unit ('ppm', 'Da', 'Th', 'm/z', etc.)
        reference_mz: Reference m/z for ppm↔Da conversion (default 500, typical peptide)

    Returns:
        Tolerance value in to_unit

    Raises:
        ValueError: If from_unit or to_unit is not recognized

    Example:
        >>> convert_tolerance(10, 'ppm', 'Da', reference_mz=500)
        0.005
        >>> convert_tolerance(0.02, 'Da', 'ppm', reference_mz=500)
        40.0
    """
    from_canon = _normalize_unit(from_unit)
    to_canon = _normalize_unit(to_unit)

    if from_canon == to_canon:
        return value
    if from_canon == "ppm" and to_canon == "Da":
        return ppm_to_da(value, reference_mz)
    if from_canon == "Da" and to_canon == "ppm":
        return da_to_ppm(value, reference_mz)
    raise ValueError(f"Cannot convert {from_unit} to {to_unit}")


def ppm_to_da(ppm: float, mz: float) -> float:
    """
    Convert ppm tolerance to Da at a given m/z.

    Args:
        ppm: Tolerance in ppm
        mz: Reference m/z value

    Returns:
        Tolerance in Da
    """
    return (ppm * mz) / 1e6


def da_to_ppm(da: float, mz: float) -> float:
    """
    Convert Da tolerance to ppm at a given m/z.

    Args:
        da: Tolerance in Da
        mz: Reference m/z value

    Returns:
        Tolerance in ppm
    """
    if mz == 0:
        return 0.0
    return (da / mz) * 1e6


def normalize_tolerance(value: float, unit: str, reference_mz: float = 500.0) -> Tuple[float, float]:
    """
    Normalize tolerance to both ppm and Da.

    Uses generic unit conversion; accepts unit aliases (ppm, Da, Th, m/z, etc.).

    Args:
        value: Tolerance value
        unit: Unit ('ppm', 'Da', 'Th', 'm/z', etc.)
        reference_mz: Reference m/z for conversion (default 500 Da, typical peptide)

    Returns:
        Tuple of (ppm_value, da_value)
    """
    ppm_val = convert_tolerance(value, unit, "ppm", reference_mz)
    da_val = convert_tolerance(value, unit, "Da", reference_mz)
    return ppm_val, da_val


# =============================================================================
# Result Data Classes
# =============================================================================

@dataclass
class ToleranceResult:
    """Result from a single tolerance estimation method."""

    method: str
    precursor_ppm: Optional[float] = None
    precursor_da: Optional[float] = None
    fragment_ppm: Optional[float] = None
    fragment_da: Optional[float] = None
    confidence: float = 0.0
    details: Dict[str, Any] = field(default_factory=dict)
    error: Optional[str] = None

    def has_precursor(self) -> bool:
        return self.precursor_ppm is not None or self.precursor_da is not None

    def has_fragment(self) -> bool:
        return self.fragment_ppm is not None or self.fragment_da is not None

    def get_precursor_ppm(self, reference_mz: float = 500.0) -> Optional[float]:
        """Get precursor tolerance in ppm, converting if necessary."""
        if self.precursor_ppm is not None:
            return self.precursor_ppm
        if self.precursor_da is not None:
            return da_to_ppm(self.precursor_da, reference_mz)
        return None

    def get_fragment_ppm(self, reference_mz: float = 500.0) -> Optional[float]:
        """Get fragment tolerance in ppm, converting if necessary."""
        if self.fragment_ppm is not None:
            return self.fragment_ppm
        if self.fragment_da is not None:
            return da_to_ppm(self.fragment_da, reference_mz)
        return None

    def get_precursor_da(self, reference_mz: float = 500.0) -> Optional[float]:
        """Get precursor tolerance in Da, converting if necessary."""
        if self.precursor_da is not None:
            return self.precursor_da
        if self.precursor_ppm is not None:
            return ppm_to_da(self.precursor_ppm, reference_mz)
        return None

    def get_fragment_da(self, reference_mz: float = 500.0) -> Optional[float]:
        """Get fragment tolerance in Da, converting if necessary."""
        if self.fragment_da is not None:
            return self.fragment_da
        if self.fragment_ppm is not None:
            return ppm_to_da(self.fragment_ppm, reference_mz)
        return None


@dataclass
class ConsensusResult:
    """Consensus result from all tolerance estimation methods."""

    # Individual method results
    gaussian: Optional[ToleranceResult] = None
    crux: Optional[ToleranceResult] = None

    # Recommended values
    recommended_precursor_ppm: Optional[float] = None
    recommended_precursor_da: Optional[float] = None
    recommended_fragment_ppm: Optional[float] = None
    recommended_fragment_da: Optional[float] = None

    # Recommendation confidence and method
    precursor_confidence: float = 0.0
    fragment_confidence: float = 0.0
    recommendation_basis: str = ""
    reference_mz: float = 500.0
    # The native unit used by the source of the fragment recommendation
    # ("ppm", "Da", or "m/z").  Useful for downstream consumers that want to
    # preserve the original unit from RunAssessor rather than converting.
    fragment_native_unit: Optional[str] = None

    def format_report(self, reference_mz: float = 500.0) -> str:
        """Format a human-readable report of all results."""
        lines = []
        lines.append("=" * 60)
        lines.append("TOLERANCE ESTIMATION REPORT")
        lines.append("=" * 60)
        lines.append(f"Reference m/z for conversions: {reference_mz} Da")
        lines.append("")

        # Precursor tolerance section
        lines.append("-" * 40)
        lines.append("PRECURSOR MASS TOLERANCE")
        lines.append("-" * 40)

        if self.gaussian and self.gaussian.has_precursor():
            ppm = self.gaussian.get_precursor_ppm(reference_mz)
            da = ppm_to_da(ppm, reference_mz) if ppm else None
            fit_quality = self.gaussian.details.get('precursor_fit_quality', 'unknown')
            lines.append(f"  Gaussian:  {ppm:.2f} ppm ({da:.4f} Da) [fit: {fit_quality}]")
        else:
            lines.append(f"  Gaussian:  not available")

        if self.crux and self.crux.has_precursor():
            ppm = self.crux.get_precursor_ppm(reference_mz)
            da = ppm_to_da(ppm, reference_mz) if ppm else None
            lines.append(f"  Crux:      {ppm:.2f} ppm ({da:.4f} Da)")
        else:
            lines.append(f"  Crux:      not available")

        if self.recommended_precursor_ppm:
            lines.append("")
            lines.append(f"  >>> RECOMMENDED: {self.recommended_precursor_ppm:.2f} ppm "
                        f"({self.recommended_precursor_da:.4f} Da)")
            lines.append(f"      Confidence: {self.precursor_confidence:.0%}")

        # Fragment tolerance section
        lines.append("")
        lines.append("-" * 40)
        lines.append("FRAGMENT MASS TOLERANCE")
        lines.append("-" * 40)

        if self.gaussian and self.gaussian.has_fragment():
            ppm = self.gaussian.get_fragment_ppm(reference_mz)
            da = ppm_to_da(ppm, reference_mz) if ppm else None
            lines.append(f"  Gaussian:  {ppm:.2f} ppm ({da:.4f} Da)")
        else:
            lines.append(f"  Gaussian:  not available")

        if self.crux and self.crux.has_fragment():
            ppm = self.crux.get_fragment_ppm(reference_mz)
            da = self.crux.fragment_da
            if da:
                lines.append(f"  Crux:      {da:.4f} Da ({ppm:.2f} ppm)")
            else:
                lines.append(f"  Crux:      {ppm:.2f} ppm")
        else:
            lines.append(f"  Crux:      not available")

        if self.recommended_fragment_ppm:
            lines.append("")
            lines.append(f"  >>> RECOMMENDED: {self.recommended_fragment_ppm:.2f} ppm "
                        f"({self.recommended_fragment_da:.4f} Da)")
            lines.append(f"      Confidence: {self.fragment_confidence:.0%}")

        lines.append("")
        lines.append("-" * 40)
        lines.append(f"Recommendation basis: {self.recommendation_basis}")
        lines.append("=" * 60)

        return "\n".join(lines)

    def format_summary(self, reference_mz: Optional[float] = None) -> str:
        """
        Format a compact summary report.

        Output format:
          Provided by Gaussian: ...
          Provided by Crux: ...
          Recommended: ...
        """
        mz_ref = reference_mz if reference_mz is not None else self.reference_mz
        lines = []

        def format_method(method: Optional[ToleranceResult]) -> str:
            if not method or (not method.has_precursor() and not method.has_fragment()):
                return "not available"
            prec_ppm = method.get_precursor_ppm(mz_ref)
            prec_da = ppm_to_da(prec_ppm, mz_ref) if prec_ppm is not None else None
            frag_ppm = method.get_fragment_ppm(mz_ref)
            frag_da = (
                method.fragment_da
                if method.fragment_da is not None
                else (ppm_to_da(frag_ppm, mz_ref) if frag_ppm is not None else None)
            )
            parts = []
            if prec_ppm is not None and prec_da is not None:
                parts.append(f"precursor {prec_ppm:.2f} ppm ({prec_da:.4f} Da)")
            if (frag_ppm is not None or frag_da is not None) and frag_da is not None:
                fppm = frag_ppm if frag_ppm is not None else da_to_ppm(frag_da, mz_ref)
                parts.append(f"fragment {fppm:.2f} ppm ({frag_da:.4f} Da)")
            return "; ".join(parts) if parts else "not available"

        lines.append(f"Provided by Gaussian: {format_method(self.gaussian)}")
        lines.append(f"Provided by Crux: {format_method(self.crux)}")

        if self.recommended_precursor_ppm is not None:
            rec_prec_da = ppm_to_da(self.recommended_precursor_ppm, mz_ref)
            rec_frag_ppm = self.recommended_fragment_ppm
            rec_frag_da = None
            if rec_frag_ppm is not None:
                rec_frag_da = ppm_to_da(rec_frag_ppm, mz_ref)
            rec_parts = [
                f"precursor {self.recommended_precursor_ppm:.2f} ppm ({rec_prec_da:.4f} Da)"
            ]
            if rec_frag_ppm is not None and rec_frag_da is not None:
                rec_parts.append(f"fragment {rec_frag_ppm:.2f} ppm ({rec_frag_da:.4f} Da)")
            lines.append(f"Recommended: {'; '.join(rec_parts)}")
        else:
            lines.append("Recommended: not available")

        return "\n".join(lines)


# =============================================================================
# Gaussian Estimator (adapted from RunAssessor)
# =============================================================================

class GaussianEstimator:
    """
    Estimates mass tolerances using RunAssessor.

    RunAssessor provides the Gaussian-based logic; local replication has been removed.
    """

    def __init__(self) -> None:
        """Initialize Gaussian estimator (RunAssessor-backed)."""

    def estimate(self, mzml_file: Path) -> ToleranceResult:
        """
        Estimate tolerances from mzML file using Gaussian fitting.

        Args:
            mzml_file: Path to mzML file

        Returns:
            ToleranceResult with estimated tolerances
        """
        result = ToleranceResult(method="gaussian")
        try:
            from runassessor import MetadataHandler, MzMLAssessor
        except ImportError:
            result.error = "RunAssessor not available"
            return result

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                metadata_path = Path(tmpdir) / "study_metadata.json"
                metadata = MetadataHandler(metadata_filepath=str(metadata_path), verbose=0)
                metadata.read_or_create()

                assessor = MzMLAssessor(str(mzml_file), metadata=metadata.metadata, verbose=0)
                assessor.read_header()
                assessor.read_spectra()
                assessor.assess_lowend_composite_spectra()
                assessor.assess_neutral_loss_composite_spectra()

                if (
                    assessor.metadata.get("state", {}).get("status") != "ERROR"
                    or "multiple fragmentation types"
                    in assessor.metadata.get("state", {}).get("message", "")
                ):
                    assessor.assess_ROIs()

                metadata.infer_search_criteria()

                file_info = metadata.metadata.get("files", {}).get(str(mzml_file), {})
                summary = file_info.get("summary", {}).get("combined summary", {})

                precursor_ppm = summary.get("recommended precursor tolerance (ppm)")
                if precursor_ppm is not None:
                    result.precursor_ppm = float(precursor_ppm)

                frag_tol = summary.get("fragmentation tolerance")
                if isinstance(frag_tol, dict):
                    frag_val = frag_tol.get("recommended fragment tolerance")
                    frag_units = frag_tol.get("recommended fragment tolerance units")
                    if frag_val is not None and frag_units:
                        if frag_units == "ppm":
                            result.fragment_ppm = float(frag_val)
                        elif frag_units in {"m/z", "Da"}:
                            result.fragment_da = float(frag_val)

                result.confidence = 0.85
                result.details = {
                    "source": "runassessor",
                    "summary": summary,
                    # Preserve RunAssessor's native fragment unit for downstream use
                    "fragment_native_unit": (
                        frag_tol.get("recommended fragment tolerance units")
                        if isinstance(frag_tol, dict) else None
                    ),
                }
                return result

        except Exception as e:
            logger.warning(f"RunAssessor estimation failed: {e}")
            result.error = str(e)
            return result


# =============================================================================
# Crux Param-Medic Estimator
# =============================================================================

class CruxEstimator:
    """
    Estimates mass tolerances using Crux param-medic tool.

    Requires crux to be installed (conda install -c bioconda crux-toolkit).
    """

    def __init__(self, charges: str = "2,3,4", min_peak_pairs: int = 20):
        """
        Initialize Crux estimator.

        Args:
            charges: Charge states to consider
            min_peak_pairs: Minimum peak pairs for estimation
        """
        self.charges = charges
        self.min_peak_pairs = min_peak_pairs

    def is_available(self) -> bool:
        """Check if crux is available on PATH."""
        import shutil
        return shutil.which("crux") is not None

    def estimate(self, mzml_file: Path) -> ToleranceResult:
        """
        Estimate tolerances using Crux param-medic.

        Args:
            mzml_file: Path to mzML file

        Returns:
            ToleranceResult with estimated tolerances
        """
        result = ToleranceResult(method="crux")

        if not self.is_available():
            result.error = "Crux not available"
            return result

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                output_dir = Path(tmpdir) / "crux-output"
                input_mzml = self._prepare_mzml(Path(mzml_file), Path(tmpdir))

                # Run param-medic
                cmd = [
                    'crux', 'param-medic',
                    str(input_mzml),
                    '--output-dir', str(output_dir),
                    '--overwrite', 'T',
                    '--pm-charges', self.charges,
                    '--pm-min-peak-pairs', str(self.min_peak_pairs),
                ]

                logger.info(f"Running: {' '.join(cmd)}")
                proc = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    encoding="utf-8",
                    errors="replace",
                    timeout=300  # 5 minute timeout
                )

                if proc.returncode != 0:
                    result.error = f"Crux failed: {proc.stderr}"
                    return result

                # Parse results from log file
                log_file = output_dir / "param-medic.log.txt"
                if log_file.exists():
                    self._parse_crux_output(log_file, result)
                else:
                    result.error = "Crux output not found"

        except subprocess.TimeoutExpired:
            result.error = "Crux timed out"
        except Exception as e:
            logger.error(f"Crux estimation failed: {e}")
            result.error = str(e)

        return result

    def _parse_crux_output(self, log_file: Path, result: ToleranceResult) -> None:
        """Parse Crux param-medic output log."""
        with open(log_file, 'r', encoding='utf-8', errors='replace') as f:
            content = f.read()

        # Parse precursor error: "INFO: Precursor error: 39.10 ppm"
        import re

        prec_match = re.search(r'Precursor error:\s*([\d.]+)\s*ppm', content)
        if prec_match:
            result.precursor_ppm = float(prec_match.group(1))
            result.confidence = 0.85

        # Parse fragment bin size: "INFO: Fragment bin size: 0.0056 Th"
        frag_match = re.search(r'Fragment bin size:\s*([\d.]+)\s*Th', content)
        if frag_match:
            result.fragment_da = float(frag_match.group(1))

        result.details['raw_output'] = content

    def _prepare_mzml(self, input_file: Path, work_dir: Path) -> Path:
        """
        Prepare mzML file for Crux param-medic.

        Creates a patched *copy* in the temp work_dir (PSI term 1003145 ->
        1000615 for Crux compatibility).  The original mzML produced by
        ThermoRawFileParser is never modified.  Uses line-by-line streaming
        so that large files don't need to be loaded entirely into memory.

        - Converts RAW -> mzML using ThermoRawFileParser.sh when needed.
        - Streams a patched copy into work_dir for Crux.
        """
        input_path = Path(input_file)
        if not input_path.exists():
            raise FileNotFoundError(str(input_path))

        if input_path.suffix.lower() == ".raw":
            mzml_path = work_dir / f"{input_path.stem}.mzML"
            cmd = [
                "ThermoRawFileParser.sh",
                f"-i={input_path}",
                f"-o={work_dir}",
                "-f=2",
            ]
            logger.info(f"Converting RAW to mzML: {' '.join(cmd)}")
            subprocess.run(
                cmd, check=True, capture_output=True, text=True,
                encoding="utf-8", errors="replace", timeout=600
            )
        else:
            mzml_path = input_path

        # Stream a patched copy into the temp dir (original stays intact).
        # Use binary mode to safely handle mzML files that may contain non-
        # UTF-8 bytes in base64-encoded spectrum data.
        fixed_mzml = work_dir / f"{mzml_path.stem}_crux.mzML"
        with open(mzml_path, "rb") as infile, \
             open(fixed_mzml, "wb") as outfile:
            for line in infile:
                outfile.write(line.replace(b"1003145", b"1000615"))

        return fixed_mzml


# =============================================================================
# Main Tolerance Estimator (Orchestrator)
# =============================================================================

class ToleranceEstimator:
    """
    Orchestrates multiple tolerance estimation methods and provides consensus.

    Combines results from:
    - Gaussian fitting (fast, no external dependencies)
    - Crux param-medic (accurate, requires crux)
    """

    def __init__(
        self,
        use_gaussian: bool = True,
        use_crux: bool = True,
        reference_mz: float = 500.0,
    ):
        """
        Initialize tolerance estimator.

        Args:
            use_gaussian: Enable Gaussian fitting method
            use_crux: Enable Crux param-medic method
            reference_mz: Reference m/z for unit conversions
        """
        self.use_gaussian = use_gaussian
        self.use_crux = use_crux
        self.reference_mz = reference_mz

        self.gaussian_estimator = GaussianEstimator() if use_gaussian else None
        self.crux_estimator = CruxEstimator() if use_crux else None

    def estimate(self, mzml_files: List[Path]) -> ConsensusResult:
        """
        Estimate tolerances using all enabled methods.

        Args:
            mzml_files: List of mzML files to analyze

        Returns:
            ConsensusResult with individual and recommended tolerances
        """
        consensus = ConsensusResult()
        consensus.reference_mz = self._compute_reference_mz(mzml_files) if mzml_files else self.reference_mz
        self.reference_mz = consensus.reference_mz

        # Run Gaussian on first file (representative)
        if self.gaussian_estimator and mzml_files:
            logger.info("Running Gaussian estimation...")
            consensus.gaussian = self.gaussian_estimator.estimate(mzml_files[0])
            if consensus.gaussian.error:
                logger.warning(f"Gaussian: {consensus.gaussian.error}")
            else:
                logger.info(f"Gaussian: precursor={consensus.gaussian.precursor_ppm} ppm")

        # Run Crux on first file
        if self.crux_estimator and mzml_files:
            logger.info("Running Crux param-medic estimation...")
            consensus.crux = self.crux_estimator.estimate(mzml_files[0])
            if consensus.crux.error:
                logger.warning(f"Crux: {consensus.crux.error}")
            else:
                logger.info(f"Crux: precursor={consensus.crux.precursor_ppm} ppm, "
                          f"fragment={consensus.crux.fragment_da} Da")

        # Compute consensus recommendations
        self._compute_consensus(consensus)

        return consensus

    def _compute_reference_mz(self, mzml_files: List[Path]) -> float:
        """
        Compute a representative precursor m/z for unit conversions.

        Uses median of observed precursor m/z values from the first file.
        Falls back to default reference_mz if pyteomics is unavailable.
        """
        try:
            from pyteomics import mzml
        except ImportError:
            logger.warning("pyteomics not available; using default reference m/z")
            return self.reference_mz

        mz_values: List[float] = []
        try:
            with mzml.read(str(mzml_files[0])) as reader:
                for spectrum in reader:
                    if "precursorList" in spectrum:
                        for precursor in spectrum.get("precursorList", {}).get("precursor", []):
                            selected_ions = precursor.get("selectedIonList", {}).get("selectedIon", [])
                            for ion in selected_ions:
                                mz = ion.get("selected ion m/z")
                                if mz:
                                    mz_values.append(float(mz))
                    if len(mz_values) >= 2000:
                        break
        except Exception as e:
            logger.warning(f"Reference m/z computation failed: {e}")
            return self.reference_mz

        if not mz_values:
            return self.reference_mz
        return float(np.median(mz_values))

    def _compute_consensus(self, consensus: ConsensusResult) -> None:
        """Compute recommended tolerances from individual results."""
        # Collect precursor estimates
        precursor_estimates: List[Tuple[float, float, str]] = []
        if consensus.gaussian and consensus.gaussian.has_precursor():
            precursor_estimates.append((
                consensus.gaussian.get_precursor_ppm(self.reference_mz),
                consensus.gaussian.confidence,
                'gaussian'
            ))
        if consensus.crux and consensus.crux.has_precursor():
            precursor_estimates.append((
                consensus.crux.get_precursor_ppm(self.reference_mz),
                consensus.crux.confidence,
                'crux'
            ))

        # Collect fragment estimates
        fragment_estimates: List[Tuple[float, float, str]] = []
        if consensus.gaussian and consensus.gaussian.has_fragment():
            fragment_estimates.append((
                consensus.gaussian.get_fragment_ppm(self.reference_mz),
                consensus.gaussian.confidence,
                'gaussian'
            ))
        if consensus.crux and consensus.crux.has_fragment():
            fragment_estimates.append((
                consensus.crux.get_fragment_ppm(self.reference_mz),
                consensus.crux.confidence,
                'crux'
            ))

        # Compute recommended precursor tolerance
        if precursor_estimates:
            rec, conf, basis = self._select_best_estimate(precursor_estimates, consensus)
            consensus.recommended_precursor_ppm = rec
            consensus.recommended_precursor_da = ppm_to_da(rec, self.reference_mz)
            consensus.precursor_confidence = conf
            consensus.recommendation_basis = basis

        # Compute recommended fragment tolerance
        if fragment_estimates:
            rec, conf, basis = self._select_best_estimate(fragment_estimates, consensus)
            consensus.recommended_fragment_ppm = rec
            # Use Crux's native Da when it provided the fragment estimate
            crux_frag_da = (
                consensus.crux.fragment_da
                if consensus.crux and consensus.crux.fragment_da is not None
                else None
            )
            consensus.recommended_fragment_da = (
                crux_frag_da
                if crux_frag_da is not None and "crux" in basis
                else ppm_to_da(rec, self.reference_mz)
            )
            consensus.fragment_confidence = conf
            if consensus.recommendation_basis:
                consensus.recommendation_basis += f"; fragment: {basis}"
            else:
                consensus.recommendation_basis = f"fragment: {basis}"

            # Track native unit from the winning method
            if "gaussian" in basis and consensus.gaussian:
                native = consensus.gaussian.details.get("fragment_native_unit")
                consensus.fragment_native_unit = native or "ppm"
            elif "crux" in basis:
                consensus.fragment_native_unit = "Da"  # Crux always reports Th/Da
            else:
                consensus.fragment_native_unit = "ppm"

    def _select_best_estimate(
        self,
        estimates: List[Tuple[float, float, str]],
        consensus: ConsensusResult
    ) -> Tuple[float, float, str]:
        """
        Select best estimate from multiple methods.

        Strategy:
        1. If all methods agree (within 20%), use weighted average
        2. If two agree and one differs, use median (robust to outlier)
        3. If all differ significantly, prefer by confidence ranking

        Args:
            estimates: List of (value, confidence, method_name) tuples

        Returns:
            Tuple of (recommended_value, confidence, basis_explanation)
        """
        if len(estimates) == 1:
            return estimates[0][0], estimates[0][1], f"single method: {estimates[0][2]}"

        values = [e[0] for e in estimates]
        confidences = [e[1] for e in estimates]
        methods = [e[2] for e in estimates]

        # Check agreement (within 20% of median)
        median_val = np.median(values)
        agreement_threshold = 0.2

        agreeing = [
            i for i, v in enumerate(values)
            if abs(v - median_val) / median_val <= agreement_threshold
        ]

        if len(agreeing) == len(values):
            # All agree - weighted average by confidence
            weights = np.array(confidences)
            weighted_avg = np.average(values, weights=weights)
            avg_confidence = np.mean(confidences)
            return (
                float(weighted_avg),
                float(avg_confidence),
                f"weighted average of {', '.join(methods)} (all agree)"
            )

        elif len(agreeing) >= 2:
            # Majority agree - use median (robust to outlier)
            return (
                float(median_val),
                float(np.mean([confidences[i] for i in agreeing])),
                f"median of agreeing methods: {[methods[i] for i in agreeing]}"
            )

        else:
            # Methods disagree significantly - fall back to highest confidence
            best_idx = int(np.argmax(confidences))
            pref_order = {"gaussian": 0, "crux": 1}
            tied = [
                i for i, c in enumerate(confidences)
                if abs(c - confidences[best_idx]) < 1e-6
            ]
            if len(tied) > 1:
                tied.sort(key=lambda i: pref_order.get(methods[i], 9))
                best_idx = tied[0]
            return (
                values[best_idx],
                confidences[best_idx],
                f"highest confidence method: {methods[best_idx]} (methods disagree)"
            )
