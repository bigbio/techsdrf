"""
PTMDetector - Detect post-translational modifications from mzML spectra.

Uses three detection tiers without database search:
  Tier 1: Reporter ion detection (TMT, iTRAQ) from MS2/MS3 spectra
  Tier 2: Diagnostic ion detection (phospho neutral losses, glyco oxonium ions)
  Tier 3: Open mass-shift detection via precursor mass pairing with
          statistical enrichment over a random null model
"""

import logging
import math
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, NamedTuple, Optional, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Maximum number of mass-shift PTMs to report per file
TOP_N_MASS_SHIFTS = 5

# Minimum observations to report a mass shift
MASS_SHIFT_MIN_OBS = 20

# Proton mass for neutral mass calculation
PROTON_MASS = 1.00728

# Quality filters for Tier 3 precursor mass collection
MIN_PEAKS_QUALITY = 15  # minimum fragment peaks for a spectrum to be used
MIN_CHARGE_QUALITY = 2  # minimum precursor charge (skip 0=unknown, 1=non-peptide)
MIN_NEUTRAL_MASS = 400.0  # Da — lower bound for peptide masses
MAX_NEUTRAL_MASS = 6000.0  # Da — upper bound for peptide masses
DEDUP_BIN_WIDTH = 0.01  # Da — bin width for unique mass deduplication


# =============================================================================
# Data classes
# =============================================================================


@dataclass
class PTMHit:
    """A single PTM detection from one file."""

    name: str
    unimod_accession: str
    tier: int  # 1=reporter, 2=diagnostic, 3=mass-shift
    evidence_count: int
    total_spectra_scanned: int
    confidence: float  # 0.0–1.0
    modification_type: str = "variable"  # "fixed" or "variable"
    mass_delta: Optional[float] = None
    enrichment: Optional[float] = None  # observed/expected ratio (tier 3)
    probability: Optional[float] = None  # 1 - p-value (tier 3)
    expected_random: Optional[float] = None  # expected count under null (tier 3)


@dataclass
class PTMDetectionResult:
    """Per-file PTM detection results."""

    mzml_file: Path
    hits: List[PTMHit] = field(default_factory=list)
    n_spectra_scanned: int = 0
    reporter_ion_type: Optional[str] = None


# =============================================================================
# Typed constants
# =============================================================================


class ReporterIonSet(NamedTuple):
    name: str
    unimod: str
    mz_values: List[float]
    min_channels: int


class DiagnosticIon(NamedTuple):
    name: str
    unimod: str
    mass_delta: float
    ion_type: str  # "neutral_loss" or "oxonium"
    target_value: float  # neutral_loss_da or oxonium_mz
    tolerance_da: float
    min_relative_intensity: float


class MassShift(NamedTuple):
    name: str
    unimod: str
    delta: float
    tolerance: float


# =============================================================================
# Constants – Tier 1: Reporter ions
# =============================================================================

_TMT6_MZ = [126.1277, 127.1311, 128.1344, 129.1378, 130.1411, 131.1382]
_TMT11_MZ = _TMT6_MZ + [127.1248, 128.1282, 129.1315, 130.1349, 131.1415]
_TMT16_MZ = _TMT11_MZ + [132.1418, 132.1449, 133.1452, 133.1483, 134.1486]
_TMT18_MZ = _TMT16_MZ + [134.1517, 135.1551]

_ITRAQ4_MZ = [114.1112, 115.1083, 116.1116, 117.1150]
_ITRAQ8_MZ = _ITRAQ4_MZ + [113.1079, 118.1120, 119.1153, 121.1220]

# Ordered largest-first so the first match is the most specific label
REPORTER_ION_SETS: List[ReporterIonSet] = [
    ReporterIonSet("TMT18plex", "UNIMOD:2016", _TMT18_MZ, 10),
    ReporterIonSet("TMT16plex", "UNIMOD:737", _TMT16_MZ, 10),
    ReporterIonSet("TMT11plex", "UNIMOD:737", _TMT11_MZ, 7),
    ReporterIonSet("TMT6plex", "UNIMOD:737", _TMT6_MZ, 4),
    ReporterIonSet("iTRAQ8plex", "UNIMOD:730", _ITRAQ8_MZ, 5),
    ReporterIonSet("iTRAQ4plex", "UNIMOD:214", _ITRAQ4_MZ, 3),
]

REPORTER_TOL_DA = 0.01

# Pre-computed unique channels for TMT plex disambiguation
_TMT6_SET = {round(m, 4) for m in _TMT6_MZ}
_TMT11_SET = {round(m, 4) for m in _TMT11_MZ}
_TMT16_SET = {round(m, 4) for m in _TMT16_MZ}

_TMT11_UNIQUE = [m for m in _TMT11_MZ if round(m, 4) not in _TMT6_SET]
_TMT16_UNIQUE = [m for m in _TMT16_MZ if round(m, 4) not in _TMT11_SET]
_TMT18_UNIQUE = [m for m in _TMT18_MZ if round(m, 4) not in _TMT16_SET]

# Mapping from label name to (unique_channels, base_channels) for disambiguation
_PLEX_DISAMBIGUATION: Dict[str, Tuple[List[float], List[float]]] = {
    "TMT18plex": (_TMT18_UNIQUE, _TMT16_MZ),
    "TMT16plex": (_TMT16_UNIQUE, _TMT11_MZ),
    "TMT11plex": (_TMT11_UNIQUE, _TMT6_MZ),
}

# Tight tolerance for unique channel verification (adjacent TMT channels
# differ by ~0.006 Da, so 0.003 Da avoids cross-matching)
_UNIQUE_CHANNEL_TOL_DA = 0.003

# Fraction of labeled spectra that must confirm unique channels
_UNIQUE_CHANNEL_MIN_FRACTION = 0.5


# =============================================================================
# Constants – Tier 2: Diagnostic ions
# =============================================================================

DIAGNOSTIC_IONS: List[DiagnosticIon] = [
    DiagnosticIon("Phospho", "UNIMOD:21", 79.9663, "neutral_loss", 97.9769, 0.02, 0.05),
    DiagnosticIon("Phospho", "UNIMOD:21", 79.9663, "neutral_loss", 79.9663, 0.02, 0.05),
    DiagnosticIon("HexNAc", "UNIMOD:43", 203.0794, "oxonium", 204.0867, 0.02, 0.05),
    DiagnosticIon("Hex(1)HexNAc(1)", "UNIMOD:43", 365.1322, "oxonium", 366.1395, 0.02, 0.02),
]


# =============================================================================
# Constants – Tier 3: Known mass shifts
# =============================================================================

KNOWN_MASS_SHIFTS: List[MassShift] = [
    MassShift("Oxidation", "UNIMOD:35", 15.9949, 0.02),
    MassShift("Deamidation", "UNIMOD:7", 0.9840, 0.02),
    MassShift("Phospho", "UNIMOD:21", 79.9663, 0.02),
    MassShift("Acetyl", "UNIMOD:1", 42.0106, 0.02),
    MassShift("Methyl", "UNIMOD:34", 14.0157, 0.02),
    MassShift("Dimethyl", "UNIMOD:36", 28.0314, 0.02),
    MassShift("Carbamidomethyl", "UNIMOD:4", 57.0215, 0.02),
    MassShift("Carbamyl", "UNIMOD:5", 43.0058, 0.02),
    MassShift("Propionamide", "UNIMOD:24", 71.0371, 0.02),
    MassShift("Pyro-glu from Q", "UNIMOD:28", -17.0265, 0.02),
    MassShift("Pyro-glu from E", "UNIMOD:27", -18.0106, 0.02),
    MassShift("Sulfo", "UNIMOD:40", 79.9568, 0.02),
    MassShift("Formyl", "UNIMOD:122", 27.9949, 0.02),
]

# UNIMOD accessions of modifications typically applied as fixed in sample prep
# (alkylation reagents). Everything else defaults to "variable".
FIXED_MODIFICATION_UNIMOD = frozenset({
    "UNIMOD:4",   # Carbamidomethyl (iodoacetamide alkylation)
    "UNIMOD:24",  # Propionamide (acrylamide alkylation)
})


# =============================================================================
# PTMDetector class
# =============================================================================


class PTMDetector:
    """Detect PTMs from mzML spectra without database search."""

    def __init__(self, max_spectra: int = 5000):
        self.max_spectra = max_spectra

    # -----------------------------------------------------------------
    # Public API
    # -----------------------------------------------------------------

    def detect(self, mzml_file: Path) -> PTMDetectionResult:
        """Run all three detection tiers on a single mzML file."""
        result = PTMDetectionResult(mzml_file=mzml_file)

        try:
            import pyopenms as oms
        except ImportError:
            logger.warning("pyopenms not available — skipping PTM detection")
            return result

        try:
            exp = oms.MSExperiment()
            oms.MzMLFile().load(str(mzml_file), exp)
        except Exception as e:
            logger.error(f"Failed to load {mzml_file} for PTM detection: {e}")
            return result

        spectra_data = self._extract_spectra(exp)
        result.n_spectra_scanned = len(spectra_data)

        if not spectra_data:
            logger.warning(f"No MS2+ spectra found in {mzml_file.name}")
            return result

        logger.info(
            f"PTM detection: scanning {len(spectra_data)} MS2+ spectra "
            f"from {mzml_file.name}..."
        )

        # Tier 1: Reporter ions
        reporter_hits = self._detect_reporter_ions(spectra_data)
        if reporter_hits:
            result.hits.extend(reporter_hits)
            result.reporter_ion_type = reporter_hits[0].name

        # Tier 2: Diagnostic ions
        result.hits.extend(self._detect_diagnostic_ions(spectra_data))

        # Tier 3: Mass shifts
        result.hits.extend(self._detect_mass_shifts(spectra_data))

        logger.info(
            f"PTM detection complete for {mzml_file.name}: "
            f"{len(result.hits)} PTM types detected"
        )
        return result

    def aggregate_results(
        self, results: List[PTMDetectionResult]
    ) -> Dict[str, Any]:
        """Aggregate PTM detection across multiple files into consensus."""
        if not results:
            return {"hits": [], "reporter_ion_type": None, "n_files": 0}

        merged: Dict[Tuple[str, int], Dict[str, Any]] = {}

        for r in results:
            for hit in r.hits:
                key = (hit.name, hit.tier)
                if key not in merged:
                    merged[key] = {
                        "name": hit.name,
                        "unimod_accession": hit.unimod_accession,
                        "tier": hit.tier,
                        "mass_delta": hit.mass_delta,
                        "modification_type": hit.modification_type,
                        "total_evidence": 0,
                        "total_scanned": 0,
                        "max_confidence": 0.0,
                        "files_detected": 0,
                        "max_enrichment": None,
                        "max_probability": None,
                        "expected_random": None,
                    }
                entry = merged[key]
                entry["total_evidence"] += hit.evidence_count
                entry["total_scanned"] += hit.total_spectra_scanned
                entry["max_confidence"] = max(entry["max_confidence"], hit.confidence)
                entry["files_detected"] += 1

                # Keep best statistical scores across files (tier 3)
                if hit.enrichment is not None:
                    prev = entry["max_enrichment"]
                    entry["max_enrichment"] = (
                        max(prev, hit.enrichment) if prev is not None else hit.enrichment
                    )
                if hit.probability is not None:
                    prev = entry["max_probability"]
                    entry["max_probability"] = (
                        max(prev, hit.probability) if prev is not None else hit.probability
                    )
                if hit.expected_random is not None:
                    prev = entry["expected_random"]
                    entry["expected_random"] = (
                        (prev + hit.expected_random) if prev is not None else hit.expected_random
                    )

        consensus_hits = sorted(
            merged.values(),
            key=lambda h: (h["tier"], -h["max_confidence"]),
        )

        reporter_counts: Counter = Counter()
        for r in results:
            if r.reporter_ion_type:
                reporter_counts[r.reporter_ion_type] += 1
        reporter_type = (
            reporter_counts.most_common(1)[0][0] if reporter_counts else None
        )

        return {
            "hits": consensus_hits,
            "reporter_ion_type": reporter_type,
            "n_files": len(results),
        }

    # -----------------------------------------------------------------
    # Spectrum extraction
    # -----------------------------------------------------------------

    def _extract_spectra(self, exp: Any) -> List[Dict[str, Any]]:
        """Extract MS2+ spectra data from a pyopenms MSExperiment."""
        spectra_data: List[Dict[str, Any]] = []
        n_scan = min(exp.getNrSpectra(), self.max_spectra)

        for i in range(n_scan):
            spec = exp.getSpectrum(i)
            if spec.getMSLevel() < 2:
                continue

            mz_array, intensity_array = spec.get_peaks()
            if len(mz_array) == 0:
                continue

            precursor_mz = None
            precursor_charge = 0
            precursors = spec.getPrecursors()
            if precursors:
                precursor_mz = precursors[0].getMZ()
                precursor_charge = precursors[0].getCharge()

            spectra_data.append({
                "mz": mz_array,
                "intensity": intensity_array,
                "precursor_mz": precursor_mz,
                "precursor_charge": precursor_charge,
            })

        return spectra_data

    # -----------------------------------------------------------------
    # Tier 1: Reporter ion detection
    # -----------------------------------------------------------------

    def _detect_reporter_ions(
        self, spectra: List[Dict[str, Any]]
    ) -> List[PTMHit]:
        """Scan low-m/z region of MS2/MS3 spectra for reporter ion clusters.

        Returns at most one PTMHit for the best-matching label type.
        """
        n_total = len(spectra)
        label_counts = self._count_reporter_matches(spectra)

        if not label_counts:
            return []

        best = self._disambiguate_plex(label_counts, n_total, spectra)
        if best is None:
            return []

        label_name, unimod_acc, count = best
        confidence = self._fraction_to_confidence(
            count / n_total,
            thresholds=[(0.5, 0.95), (0.2, 0.85), (0.1, 0.7), (0.0, 0.5)],
        )

        return [PTMHit(
            name=label_name,
            unimod_accession=unimod_acc,
            tier=1,
            evidence_count=count,
            total_spectra_scanned=n_total,
            confidence=confidence,
            modification_type="fixed",
        )]

    @staticmethod
    def _count_reporter_matches(
        spectra: List[Dict[str, Any]],
    ) -> Dict[str, int]:
        """Count spectra matching each reporter ion set."""
        counts: Dict[str, int] = {}

        for sp in spectra:
            mz = sp["mz"]
            intensity = sp["intensity"]

            mask = (mz >= 100.0) & (mz <= 140.0)
            low_mz = mz[mask]
            low_int = intensity[mask]

            if len(low_mz) == 0:
                continue

            for rset in REPORTER_ION_SETS:
                matched = 0
                for ref in rset.mz_values:
                    diffs = np.abs(low_mz - ref)
                    if np.any(diffs <= REPORTER_TOL_DA):
                        idx = np.argmin(diffs)
                        if low_int[idx] > 0:
                            matched += 1
                if matched >= rset.min_channels:
                    counts[rset.name] = counts.get(rset.name, 0) + 1

        return counts

    def _disambiguate_plex(
        self,
        label_counts: Dict[str, int],
        n_total: int,
        spectra: List[Dict[str, Any]],
    ) -> Optional[Tuple[str, str, int]]:
        """Select the correct TMT plex level using unique channel validation.

        Iterates from most specific to least specific. For each candidate,
        checks that channels unique to that plex are truly present (not just
        noise or isotope peaks from a lower plex).
        """
        for rset in REPORTER_ION_SETS:
            base_count = label_counts.get(rset.name, 0)
            if base_count == 0 or base_count / n_total < 0.05:
                continue

            disambiguation = _PLEX_DISAMBIGUATION.get(rset.name)
            if disambiguation is not None:
                unique_mzs, base_mzs = disambiguation
                confirmed = self._count_unique_channel_spectra(
                    spectra, unique_mzs, base_mzs
                )
                if confirmed < base_count * _UNIQUE_CHANNEL_MIN_FRACTION:
                    continue

            return (rset.name, rset.unimod, base_count)

        return None

    @staticmethod
    def _count_unique_channel_spectra(
        spectra: List[Dict[str, Any]],
        unique_mzs: List[float],
        base_mzs: List[float],
    ) -> int:
        """Count spectra where unique channels are confirmed.

        Uses tight tolerance to avoid confusing TMT plex levels when
        adjacent channels differ by only ~0.006 Da. Verifies that matched
        peaks are closer to the unique channel than to any base-set channel.
        """
        if not unique_mzs:
            return 0

        min_required = max(2, (len(unique_mzs) + 1) // 2)
        base_arr = np.array(base_mzs) if base_mzs else np.array([])
        count = 0

        for sp in spectra:
            mz = sp["mz"]
            mask = (mz >= 100.0) & (mz <= 140.0)
            low_mz = mz[mask]
            if len(low_mz) == 0:
                continue

            n_matched = 0
            for ref in unique_mzs:
                diffs = np.abs(low_mz - ref)
                best_idx = np.argmin(diffs)
                if diffs[best_idx] > _UNIQUE_CHANNEL_TOL_DA:
                    continue
                if len(base_arr) > 0:
                    base_diffs = np.abs(low_mz[best_idx] - base_arr)
                    if np.min(base_diffs) <= diffs[best_idx]:
                        continue
                n_matched += 1

            if n_matched >= min_required:
                count += 1

        return count

    # -----------------------------------------------------------------
    # Tier 2: Diagnostic ions (neutral losses & oxonium ions)
    # -----------------------------------------------------------------

    def _detect_diagnostic_ions(
        self, spectra: List[Dict[str, Any]]
    ) -> List[PTMHit]:
        """Look for neutral losses and oxonium ions indicative of PTMs."""
        n_total = len(spectra)
        diag_counts: Dict[str, int] = {}
        diag_meta: Dict[str, DiagnosticIon] = {}

        for sp in spectra:
            mz = sp["mz"]
            intensity = sp["intensity"]
            precursor_mz = sp["precursor_mz"]
            charge = sp["precursor_charge"] or 2

            if len(mz) == 0:
                continue

            max_int = np.max(intensity)
            if max_int == 0:
                continue

            already_counted: set = set()

            for diag in DIAGNOSTIC_IONS:
                if diag.name in already_counted:
                    continue

                target = self._compute_diagnostic_target(diag, precursor_mz, charge)
                if target is None:
                    continue

                diffs = np.abs(mz - target)
                best_idx = np.argmin(diffs)
                if (
                    diffs[best_idx] <= diag.tolerance_da
                    and intensity[best_idx] / max_int >= diag.min_relative_intensity
                ):
                    diag_counts[diag.name] = diag_counts.get(diag.name, 0) + 1
                    diag_meta[diag.name] = diag
                    already_counted.add(diag.name)

        hits = []
        for name, count in diag_counts.items():
            fraction = count / n_total
            if fraction < 0.01:
                continue

            meta = diag_meta[name]
            confidence = self._fraction_to_confidence(
                fraction,
                thresholds=[(0.2, 0.9), (0.1, 0.7), (0.05, 0.5), (0.0, 0.3)],
            )
            hits.append(PTMHit(
                name=meta.name,
                unimod_accession=meta.unimod,
                tier=2,
                evidence_count=count,
                total_spectra_scanned=n_total,
                confidence=confidence,
                mass_delta=meta.mass_delta,
            ))

        return hits

    @staticmethod
    def _compute_diagnostic_target(
        diag: DiagnosticIon,
        precursor_mz: Optional[float],
        charge: int,
    ) -> Optional[float]:
        """Compute the expected m/z for a diagnostic ion."""
        if diag.ion_type == "neutral_loss":
            if precursor_mz is None:
                return None
            target = precursor_mz - diag.target_value / charge
            return target if target > 0 else None
        return diag.target_value

    # -----------------------------------------------------------------
    # Tier 3: Open mass-shift detection
    # -----------------------------------------------------------------

    def _detect_mass_shifts(
        self, spectra: List[Dict[str, Any]]
    ) -> List[PTMHit]:
        """Detect PTMs via precursor mass pairing with statistical scoring.

        For each known PTM mass shift, counts precursor mass pairs
        separated by that shift, then scores against a uniform-density
        null model using Poisson statistics.

        Quality filtering and deduplication are applied to reduce noise:
        - Spectra with < MIN_PEAKS_QUALITY fragment peaks are excluded
        - Spectra with charge < MIN_CHARGE_QUALITY are excluded
        - Neutral masses outside MIN_NEUTRAL_MASS–MAX_NEUTRAL_MASS are excluded
        - Masses are deduplicated to unique species (binned at DEDUP_BIN_WIDTH Da)
        """
        all_masses, unique_masses = self._collect_precursor_masses(spectra)

        logger.debug(
            f"Tier 3 quality filter: {len(spectra)} spectra -> "
            f"{len(all_masses)} passed -> {len(unique_masses)} unique masses"
        )

        if len(unique_masses) < 10:
            return []

        n_unique = len(unique_masses)
        masses = np.sort(unique_masses)

        mass_range = masses[-1] - masses[0]
        if mass_range <= 0:
            return []
        density = n_unique / mass_range

        candidates = []
        for shift in KNOWN_MASS_SHIFTS:
            count = self._count_mass_pairs(masses, abs(shift.delta), shift.tolerance)
            expected = max(n_unique * density * 2 * shift.tolerance, 0.1)
            enrichment = count / expected
            p_value = self._poisson_sf(count, expected) if count > 0 else 1.0

            candidates.append({
                "shift": shift,
                "count": count,
                "expected": expected,
                "enrichment": enrichment,
                "p_value": p_value,
                "probability": 1.0 - p_value,
            })

        significant = [
            c for c in candidates
            if c["enrichment"] > 1.5 and c["p_value"] < 0.01 and c["count"] >= MASS_SHIFT_MIN_OBS
        ]
        significant.sort(key=lambda c: -c["enrichment"])

        hits = []
        for c in significant[:TOP_N_MASS_SHIFTS]:
            shift = c["shift"]
            enrichment = c["enrichment"]
            confidence = self._fraction_to_confidence(
                enrichment,
                thresholds=[(5.0, 0.95), (3.0, 0.85), (2.0, 0.7), (0.0, 0.5)],
            )
            mod_type = (
                "fixed" if shift.unimod in FIXED_MODIFICATION_UNIMOD else "variable"
            )
            hits.append(PTMHit(
                name=shift.name,
                unimod_accession=shift.unimod,
                tier=3,
                evidence_count=c["count"],
                total_spectra_scanned=n_unique,
                confidence=confidence,
                modification_type=mod_type,
                mass_delta=shift.delta,
                enrichment=round(enrichment, 2),
                probability=round(c["probability"], 4),
                expected_random=round(c["expected"], 1),
            ))

        return hits

    @staticmethod
    def _collect_precursor_masses(
        spectra: List[Dict[str, Any]],
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Collect and filter precursor neutral masses from MS2 spectra.

        Applies quality gates:
        - Minimum fragment peak count (MIN_PEAKS_QUALITY)
        - Minimum charge state (MIN_CHARGE_QUALITY)
        - Neutral mass within peptide range (MIN/MAX_NEUTRAL_MASS)

        Then deduplicates by binning to DEDUP_BIN_WIDTH Da resolution.

        Returns:
            Tuple of (all_filtered_masses, unique_masses) as numpy arrays.
        """
        masses = []
        for sp in spectra:
            mz = sp["precursor_mz"]
            if mz is None or mz <= 0:
                continue
            charge = sp["precursor_charge"]
            # Skip unknown charge (0) and singly-charged (non-peptide in ESI)
            if charge < MIN_CHARGE_QUALITY:
                continue
            # Skip spectra with too few fragment peaks
            if len(sp["mz"]) < MIN_PEAKS_QUALITY:
                continue
            neutral_mass = (mz - PROTON_MASS) * charge
            # Skip masses outside peptide range
            if neutral_mass < MIN_NEUTRAL_MASS or neutral_mass > MAX_NEUTRAL_MASS:
                continue
            masses.append(neutral_mass)

        all_masses = np.array(masses) if masses else np.array([], dtype=float)

        # Deduplicate: bin to DEDUP_BIN_WIDTH resolution, keep unique
        if len(all_masses) > 0:
            binned = np.round(all_masses / DEDUP_BIN_WIDTH) * DEDUP_BIN_WIDTH
            unique_masses = np.unique(binned)
        else:
            unique_masses = np.array([], dtype=float)

        return all_masses, unique_masses

    @staticmethod
    def _count_mass_pairs(
        sorted_masses: np.ndarray, delta: float, tolerance: float
    ) -> int:
        """Count masses M where M+delta exists in the sorted array."""
        count = 0
        targets = sorted_masses + delta
        for target in targets:
            idx = np.searchsorted(sorted_masses, target - tolerance)
            if idx < len(sorted_masses) and abs(sorted_masses[idx] - target) <= tolerance:
                count += 1
            elif idx > 0 and abs(sorted_masses[idx - 1] - target) <= tolerance:
                count += 1
        return count

    # -----------------------------------------------------------------
    # Utilities
    # -----------------------------------------------------------------

    @staticmethod
    def _fraction_to_confidence(
        value: float, thresholds: List[Tuple[float, float]]
    ) -> float:
        """Map a value to a confidence score using threshold brackets.

        Thresholds are (min_value, confidence) pairs, checked in order.
        Returns the confidence for the first matching threshold.
        """
        for threshold, confidence in thresholds:
            if value >= threshold:
                return confidence
        return thresholds[-1][1]

    @staticmethod
    def _poisson_sf(k: int, lam: float) -> float:
        """Poisson survival function P(X >= k).

        Uses scipy when available, with fallback to normal approximation
        (large lambda) or direct log-space summation (small lambda).
        """
        try:
            from scipy.stats import poisson
            return float(poisson.sf(k - 1, lam))
        except ImportError:
            pass

        if lam > 100:
            z = (k - lam) / math.sqrt(lam)
            return 0.5 * math.erfc(z / math.sqrt(2))

        # Direct summation: P(X >= k) = 1 - sum_{i=0}^{k-1} e^{-lam} * lam^i / i!
        cum = 0.0
        log_lam = math.log(lam) if lam > 0 else -float("inf")
        for i in range(k):
            log_prob = -lam + i * log_lam - math.lgamma(i + 1)
            cum += math.exp(log_prob)
        return max(0.0, 1.0 - cum)
