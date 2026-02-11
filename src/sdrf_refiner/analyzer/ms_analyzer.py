"""
MSAnalyzer - Analyze mzML files to detect instrument parameters.

Adapts algorithms from RunAssessor for instrument detection, fragmentation
type parsing, and tolerance calculation.
"""

import gzip
import logging
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
from sdrf_refiner.config import (
    ACQUISITION_ONTOLOGY,
    FRAGMENTATION_ONTOLOGY,
    IONIZATION_ONTOLOGY,
    MASS_ANALYZER_ONTOLOGY,
    format_ontology_value,
)

logger = logging.getLogger(__name__)


@dataclass
class AnalysisResult:
    """Results from analyzing an mzML file."""

    mzml_file: Path
    instrument_model: Dict[str, str] = field(default_factory=dict)
    instrument_category: str = "unknown"
    fragmentation_type: str = "unknown"  # Primary/MS2 fragmentation
    ms2_fragmentation_type: str = "unknown"  # MS2 specific
    ms3_fragmentation_type: str = "unknown"  # MS3 specific (for SPS-MS3, etc.)
    high_accuracy_precursors: bool = False
    precursor_tolerance_ppm: Optional[float] = None
    fragment_tolerance_ppm: Optional[float] = None
    fragment_tolerance_da: Optional[float] = None
    acquisition_type: str = "unknown"
    n_spectra: int = 0
    n_ms1_spectra: int = 0
    n_ms2_spectra: int = 0
    n_ms3_spectra: int = 0  # MS3 spectrum count
    ms_min_charge: Optional[int] = None  # Minimum observed precursor charge
    ms_max_charge: Optional[int] = None  # Maximum observed precursor charge
    collision_energy: Optional[str] = None  # e.g. "32 eV" or "30 NCE" (SDRF format)
    # --- Properties detected via pyopenms header parsing ---
    mass_analyzer_type: str = "unknown"  # e.g. "Orbitrap", "TOF", "IT", "Quadrupole"
    ionization_method: str = "unknown"  # e.g. "ESI", "NESI", "MALDI"
    polarity: str = "unknown"  # "positive", "negative", "mixed"
    ms1_scan_range: Optional[tuple] = None  # (low_mz, high_mz)
    isolation_window_da: Optional[float] = None  # Median full isolation width in Da
    instrument_serial_number: Optional[str] = None  # Instrument serial if available
    # Tolerance source: "heuristic" (instrument category lookup) or "empirical" (Crux/Gaussian)
    precursor_tolerance_source: str = "heuristic"
    fragment_tolerance_source: str = "heuristic"
    raw_stats: Dict[str, Any] = field(default_factory=dict)
    confidence: Dict[str, float] = field(default_factory=dict)


class MSAnalyzer:
    """
    Analyzes mzML files to detect MS instrument parameters.

    Uses adapted algorithms from RunAssessor for:
    - Instrument detection via CV terms in mzML header
    - Fragmentation type detection via Thermo filter strings
    - Precursor/fragment tolerance estimation
    """

    # Instrument categories by CV accession
    INSTRUMENT_BY_CATEGORY = {
        "pureHCD": [
            ("MS:1000649", "Exactive"),
            ("MS:1002526", "Exactive Plus"),
            ("MS:1001911", "Q Exactive"),
            ("MS:1002993", "Q Exactive Focus"),
            ("MS:1002523", "Q Exactive HF"),
            ("MS:1002877", "Q Exactive HF-X"),
            ("MS:1002634", "Q Exactive Plus"),
            ("MS:1003245", "Q Exactive UHMR"),
            ("MS:1003378", "Orbitrap Astral"),
        ],
        "ion_trap": [
            ("MS:1000447", "LTQ"),
            ("MS:1000855", "LTQ Velos"),
            ("MS:1000854", "LTQ XL"),
            ("MS:1000167", "LCQ Advantage"),
            ("MS:1000168", "LCQ Classic"),
            ("MS:1001909", "Velos Plus"),
        ],
        "variable": [
            ("MS:1000448", "LTQ FT"),
            ("MS:1000557", "LTQ FT Ultra"),
            ("MS:1000449", "LTQ Orbitrap"),
            ("MS:1000555", "LTQ Orbitrap Discovery"),
            ("MS:1001742", "LTQ Orbitrap Velos"),
            ("MS:1000556", "LTQ Orbitrap XL"),
            ("MS:1003029", "Orbitrap Eclipse"),
            ("MS:1001910", "Orbitrap Elite"),
            ("MS:1003095", "Orbitrap Exploris 120"),
            ("MS:1003094", "Orbitrap Exploris 240"),
            ("MS:1003028", "Orbitrap Exploris 480"),
            ("MS:1002416", "Orbitrap Fusion"),
            ("MS:1002732", "Orbitrap Fusion Lumos"),
        ],
        "QTOF": [
            ("MS:1000126", "Waters instrument model"),
            ("MS:1000122", "Bruker Daltonics instrument model"),
            ("MS:1001547", "Bruker Daltonics maXis series"),
            ("MS:1003123", "Bruker Daltonics timsTOF series"),
            ("MS:1003005", "timsTOF Pro"),
            ("MS:1003230", "timsTOF Pro 2"),
            ("MS:1000121", "AB SCIEX instrument model"),
            ("MS:1002583", "TripleTOF 4600"),
            ("MS:1000932", "TripleTOF 5600"),
            ("MS:1002584", "TripleTOF 5600+"),
            ("MS:1002533", "TripleTOF 6600"),
            ("MS:1003293", "ZenoTOF 7600"),
        ],
    }

    def __init__(self, verbose: int = 0):
        """
        Initialize MS analyzer.

        Args:
            verbose: Verbosity level (0=quiet, 1=info, 2=debug)
        """
        self.verbose = verbose
        self._pyopenms_extras: Dict[str, Any] = {}  # Filled by _read_header_pyopenms

        # Build instrument lookup dict
        self._instrument_lookup = {}
        for category, instruments in self.INSTRUMENT_BY_CATEGORY.items():
            for accession, name in instruments:
                self._instrument_lookup[accession] = {
                    "name": name,
                    "accession": accession,
                    "category": category,
                }

    def analyze(self, mzml_file: Path) -> AnalysisResult:
        """
        Analyze an mzML file and extract instrument parameters.

        Args:
            mzml_file: Path to mzML file

        Returns:
            AnalysisResult with detected parameters
        """
        mzml_file = Path(mzml_file)
        logger.info(f"Analyzing {mzml_file.name}...")

        result = AnalysisResult(mzml_file=mzml_file)

        try:
            # Step 1: Read header for instrument detection
            self._pyopenms_extras = {}  # Reset per file
            instrument_data = self._read_header(mzml_file)
            if instrument_data:
                result.instrument_model = instrument_data
                result.instrument_category = instrument_data.get("category", "unknown")
                result.confidence["instrument"] = 0.95

            # Merge pyopenms extras (mass analyzer, ionization, polarity, etc.)
            extras = self._pyopenms_extras
            if extras:
                result.mass_analyzer_type = extras.get("mass_analyzer_type", "unknown")
                result.ionization_method = extras.get("ionization_method", "unknown")
                result.polarity = extras.get("polarity", "unknown")
                result.ms1_scan_range = extras.get("ms1_scan_range")
                result.isolation_window_da = extras.get("isolation_window_da")
                result.instrument_serial_number = extras.get("instrument_serial_number")

                # Use mass analyzer to refine instrument category if needed
                if result.instrument_category == "unknown":
                    ma = result.mass_analyzer_type
                    if ma in ("Orbitrap", "FourierTransform"):
                        result.instrument_category = "Orbitrap"
                    elif ma == "TOF":
                        result.instrument_category = "QTOF"
                    elif ma in ("IT", "LIT", "PaulIonTrap"):
                        result.instrument_category = "IT"

            # Step 2: Read spectra for fragmentation type and stats
            stats = self._read_spectra(mzml_file)
            result.fragmentation_type = stats.get("fragmentation_type", "unknown")
            result.ms2_fragmentation_type = stats.get("ms2_fragmentation_type", "unknown")
            result.ms3_fragmentation_type = stats.get("ms3_fragmentation_type", "unknown")
            result.high_accuracy_precursors = stats.get("high_accuracy_precursors") == "true"
            result.n_spectra = stats.get("n_spectra", 0)
            result.n_ms1_spectra = stats.get("n_ms1_spectra", 0)
            result.n_ms2_spectra = stats.get("n_ms2_spectra", 0)
            result.n_ms3_spectra = stats.get("n_ms3_spectra", 0)
            result.raw_stats = stats

            # Extract collision energy - spec: "30 NCE", "30% NCE", "25 NCE;27 NCE;30 NCE"
            ce_dict = stats.get("collision_energies", {})
            if ce_dict:
                total = sum(ce_dict.values())
                # Stepped CE: 2+ distinct values with >5% count each
                candidates = [
                    (k, c) for k, c in ce_dict.items() if c >= max(2, total * 0.05)
                ]
                if len(candidates) >= 2:
                    # Output stepped format: "25 NCE;27 NCE;30 NCE" (sorted by value)
                    def parse_ce(s: str) -> float:
                        try:
                            return float(s.split()[0].rstrip("%"))
                        except (ValueError, IndexError, AttributeError):
                            return 0.0

                    sorted_ce = sorted(
                        [k for k, _ in candidates],
                        key=parse_ce,
                    )
                    result.collision_energy = ";".join(sorted_ce)
                    result.confidence["collision_energy"] = 0.85
                else:
                    dominant_ce = max(ce_dict, key=ce_dict.get)
                    result.collision_energy = dominant_ce
                    result.confidence["collision_energy"] = 0.9

            # Extract min/max charge states from stats
            charge_states = []
            for key in stats:
                if key.startswith("n_charge_") and key.endswith("_precursors"):
                    # Extract charge value from key like "n_charge_2_precursors"
                    try:
                        charge = int(key.split("_")[2])
                        if stats[key] > 0:  # Only count if we have precursors with this charge
                            charge_states.append(charge)
                    except (ValueError, IndexError):
                        pass
            if charge_states:
                result.ms_min_charge = min(charge_states)
                result.ms_max_charge = max(charge_states)

            # Determine acquisition type from isolation windows
            if "isolation_window_full_widths" in stats:
                widths = stats["isolation_window_full_widths"]
                if widths:
                    median_width = np.median(list(widths.keys()))
                    if median_width <= 3.0:
                        result.acquisition_type = "DDA"
                    elif median_width >= 15.0:
                        result.acquisition_type = "DIA"
                    else:
                        result.acquisition_type = "unknown"

            # Step 3: Estimate tolerances based on instrument category
            self._estimate_tolerances(result)

            if result.fragmentation_type != "unknown":
                result.confidence["fragmentation"] = 0.9

            logger.info(
                f"Analysis complete: instrument={result.instrument_model.get('name', 'unknown')}, "
                f"fragmentation={result.fragmentation_type}, "
                f"spectra={result.n_spectra}"
            )

        except Exception as e:
            logger.error(f"Analysis failed for {mzml_file}: {e}")
            result.raw_stats["error"] = str(e)

        return result

    def _read_header(self, mzml_file: Path) -> Optional[Dict[str, str]]:
        """
        Read mzML header to detect instrument model via pyopenms.

        Args:
            mzml_file: Path to mzML file

        Returns:
            Instrument model data or None
        """
        return self._read_header_pyopenms(mzml_file)

    # --- pyopenms enum → human-readable mappings ---
    _MASS_ANALYZER_NAMES = {
        0: "unknown", 1: "Quadrupole", 2: "PaulIonTrap", 3: "RadialEjectionLIT",
        4: "AxialEjectionLIT", 5: "TOF", 6: "Sector", 7: "FourierTransform",
        8: "IonStorage", 9: "ESA", 10: "IT", 11: "SWIFT", 12: "Cyclotron",
        13: "Orbitrap", 14: "LIT",
    }
    _IONIZATION_NAMES = {
        0: "unknown", 1: "ESI", 2: "EI", 3: "CI", 4: "FAB", 5: "TSP",
        6: "LD", 7: "FD", 8: "FI", 9: "PD", 10: "SI", 11: "TI", 12: "API",
        13: "ISI", 14: "CID", 15: "CAD", 16: "HN", 17: "APCI", 18: "APPI",
        19: "ICP", 20: "NESI", 21: "MESI", 22: "SELDI", 23: "SEND", 24: "FIB",
        25: "MALDI", 26: "MPI", 27: "DI", 28: "FA", 29: "FII", 30: "GD_MS",
        31: "NICI", 32: "NRMS", 33: "PI", 34: "PYMS", 35: "REMPI", 36: "AI",
        37: "ASI", 38: "AD", 39: "AUI", 40: "CEI", 41: "CHEMI", 42: "DISSI",
        43: "LSI", 44: "PEI", 45: "SOI", 46: "SPI", 47: "SUI", 48: "VI",
        49: "AP_MALDI", 50: "SILI", 51: "SALDI",
    }
    _POLARITY_NAMES = {0: "unknown", 1: "positive", 2: "negative"}

    def _read_header_pyopenms(self, mzml_file: Path) -> Optional[Dict[str, str]]:
        """
        Read instrument model and additional metadata from mzML using pyopenms.

        In addition to instrument model, this method populates:
        - mass_analyzer_type, ionization_method, polarity
        - ms1_scan_range, isolation_window_da, instrument_serial_number

        These are stored in ``self._pyopenms_extras`` so the caller
        (``analyze``) can merge them into the ``AnalysisResult``.

        Returns:
            Instrument model data, or None if pyopenms is unavailable / fails
        """
        try:
            import pyopenms as oms
        except ImportError:
            logger.error("pyopenms is required but not installed")
            return None

        try:
            exp = oms.MSExperiment()
            oms.MzMLFile().load(str(mzml_file), exp)

            instrument = exp.getInstrument()
            inst_name = str(instrument.getName()).strip()
            inst_model = str(instrument.getModel()).strip()

            # ---- Extra properties from instrument configuration ----
            extras: Dict[str, Any] = {}

            # Mass analyzer type
            mass_analyzers = instrument.getMassAnalyzers()
            if mass_analyzers:
                ma_type_int = mass_analyzers[0].getType()
                extras["mass_analyzer_type"] = self._MASS_ANALYZER_NAMES.get(
                    ma_type_int, f"type_{ma_type_int}"
                )

            # Ion source / ionization method
            ion_sources = instrument.getIonSources()
            if ion_sources:
                ionization_int = ion_sources[0].getIonizationMethod()
                extras["ionization_method"] = self._IONIZATION_NAMES.get(
                    ionization_int, f"method_{ionization_int}"
                )

            # Instrument serial number (metavalue)
            keys: List = []
            instrument.getKeys(keys)
            for k in keys:
                k_str = k.decode() if isinstance(k, bytes) else str(k)
                if "serial" in k_str.lower():
                    extras["instrument_serial_number"] = str(
                        instrument.getMetaValue(k)
                    )
                    break

            # ---- Spectrum-level aggregated properties ----
            polarity_counts: Dict[int, int] = {}
            scan_windows: List[tuple] = []
            isolation_widths: List[float] = []

            n_sample = min(exp.getNrSpectra(), 5000)
            for i in range(n_sample):
                spec = exp.getSpectrum(i)
                ml = spec.getMSLevel()
                isettings = spec.getInstrumentSettings()

                # Polarity
                pol = isettings.getPolarity()
                polarity_counts[pol] = polarity_counts.get(pol, 0) + 1

                # MS1 scan window
                if ml == 1 and not scan_windows:
                    for w in isettings.getScanWindows():
                        scan_windows.append((w.begin, w.end))

                # Isolation window from precursors
                if ml >= 2:
                    for prec in spec.getPrecursors():
                        lo = prec.getIsolationWindowLowerOffset()
                        up = prec.getIsolationWindowUpperOffset()
                        if lo > 0 or up > 0:
                            isolation_widths.append(lo + up)

            # Polarity consensus
            if polarity_counts:
                # Remove "unknown" (0) if we have real polarity data
                known = {k: v for k, v in polarity_counts.items() if k != 0}
                if known:
                    if len(known) == 1:
                        extras["polarity"] = self._POLARITY_NAMES.get(
                            next(iter(known)), "unknown"
                        )
                    else:
                        extras["polarity"] = "mixed"

            # MS1 scan range
            if scan_windows:
                extras["ms1_scan_range"] = scan_windows[0]

            # Median isolation window
            if isolation_widths:
                import statistics
                extras["isolation_window_da"] = round(
                    statistics.median(isolation_widths), 2
                )

            # Store extras for the caller
            self._pyopenms_extras = extras

            # ---- Instrument model matching ----
            candidate_names = [inst_name, inst_model]

            # Build a reverse lookup: instrument name -> accession data
            name_lookup = {}
            for accession, data in self._instrument_lookup.items():
                name_lookup[data["name"].lower()] = data

            for name in candidate_names:
                if not name:
                    continue
                key = name.lower()
                if key in name_lookup:
                    model_data = name_lookup[key].copy()
                    logger.info(f"pyopenms: detected instrument {model_data['name']} "
                                f"({model_data['accession']})")
                    return model_data

                # Partial match: check if any known name is contained
                for known_name, data in name_lookup.items():
                    if known_name in key or key in known_name:
                        model_data = data.copy()
                        logger.info(f"pyopenms: detected instrument {model_data['name']} "
                                    f"({model_data['accession']}) via partial match")
                        return model_data

            # pyopenms loaded the file but we couldn't match to a known instrument
            display_name = inst_name or inst_model or "unknown"
            logger.warning(f"pyopenms: instrument '{display_name}' not in known instrument list")
            # Use mass analyzer type to infer category
            ma = extras.get("mass_analyzer_type", "unknown")
            if ma in ("Orbitrap", "FourierTransform"):
                category = "Orbitrap"
            elif ma == "TOF":
                category = "QTOF"
            elif ma in ("IT", "LIT", "PaulIonTrap", "RadialEjectionLIT", "AxialEjectionLIT"):
                category = "IT"
            else:
                category = "unknown"
            return {"name": display_name, "accession": "", "category": category}

        except Exception as e:
            logger.warning(f"pyopenms header parsing failed: {e}")
            return None

    def _read_spectra(self, mzml_file: Path) -> Dict[str, Any]:
        """
        Read spectra to determine fragmentation type and collect stats.

        Args:
            mzml_file: Path to mzML file

        Returns:
            Statistics dictionary
        """
        try:
            from pyteomics import mzml
        except ImportError:
            logger.warning("pyteomics not available, using limited spectrum analysis")
            return self._read_spectra_basic(mzml_file)

        stats = {
            "n_spectra": 0,
            "n_ms1_spectra": 0,
            "n_ms2_spectra": 0,
            "n_ms3_spectra": 0,
            "fragmentation_type": "unknown",
            "ms2_fragmentation_type": "unknown",
            "ms3_fragmentation_type": "unknown",
            "high_accuracy_precursors": "unknown",
            "isolation_window_full_widths": {},
            "collision_energies": {},  # value+unit -> count
        }

        # Initialize fragmentation counters for MS2 and MS3
        for ms_level in ["ms2", "ms3"]:
            for prefix in ["HR", "LR"]:
                for frag in ["HCD", "IT_CID", "IT_ETD", "EThcD", "ETciD"]:
                    stats[f"n_{ms_level}_{prefix}_{frag}_spectra"] = 0

        try:
            with mzml.read(str(mzml_file)) as reader:
                for spectrum in reader:
                    stats["n_spectra"] += 1

                    # Check for filter string (Thermo data)
                    filter_string = spectrum.get("filter string")
                    if filter_string:
                        self._parse_filter_string(filter_string, stats)
                    else:
                        # Non-Thermo data - use ms level
                        ms_level = spectrum.get("ms level", 1)
                        if ms_level == 1:
                            stats["n_ms1_spectra"] += 1
                        elif ms_level == 2:
                            stats["n_ms2_spectra"] += 1
                            # Try to determine fragmentation from CV terms
                            self._detect_fragmentation_cv(spectrum, stats)

                    # Record collision energy from activation (MS:1000045)
                    if "precursorList" in spectrum:
                        for precursor in spectrum.get("precursorList", {}).get(
                            "precursor", []
                        ):
                            act = precursor.get("activation", {})
                            # MS:1000045 collision energy; MS:1000138 percent collision energy
                            ce_val = act.get("collision energy")
                            pct_val = act.get("percent collision energy")
                            if pct_val is not None:
                                try:
                                    v = float(pct_val)
                                    key = f"{int(v) if v == int(v) else v}% NCE"
                                    stats["collision_energies"][key] = (
                                        stats["collision_energies"].get(key, 0) + 1
                                    )
                                except (TypeError, ValueError):
                                    pass
                            elif ce_val is not None:
                                # MS:1000045 - mzML uses eV; Thermo often reports NCE
                                # Spec: "30 NCE", "27 eV", "30% NCE"
                                try:
                                    v = float(ce_val)
                                    unit = "NCE" if 15 <= v <= 45 else "eV"
                                    key = f"{int(v) if v == int(v) else v} {unit}"
                                    stats["collision_energies"][key] = (
                                        stats["collision_energies"].get(key, 0) + 1
                                    )
                                except (TypeError, ValueError):
                                    pass
                            # Isolation window
                            isolation = precursor.get("isolationWindow", {})
                            lower = isolation.get(
                                "isolation window lower offset", 0
                            )
                            upper = isolation.get(
                                "isolation window upper offset", 0
                            )
                            if lower and upper:
                                width = float(lower) + float(upper)
                                stats["isolation_window_full_widths"][width] = (
                                    stats["isolation_window_full_widths"].get(width, 0)
                                    + 1
                                )

                            # Charge state from selected ion
                            selected_ions = precursor.get("selectedIonList", {}).get(
                                "selectedIon", []
                            )
                            for ion in selected_ions:
                                charge = ion.get("charge state")
                                if charge is not None:
                                    try:
                                        charge = int(charge)
                                        charge_key = f"n_charge_{charge}_precursors"
                                        stats[charge_key] = stats.get(charge_key, 0) + 1
                                    except (ValueError, TypeError):
                                        pass

                    # Limit number of spectra for performance
                    if stats["n_spectra"] >= 10000:
                        break

        except Exception as e:
            logger.warning(f"Error reading spectra: {e}")

        # Determine dominant fragmentation type
        self._determine_fragmentation_type(stats)

        return stats

    def _read_spectra_basic(self, mzml_file: Path) -> Dict[str, Any]:
        """Basic spectrum reading without pyteomics."""
        stats = {
            "n_spectra": 0,
            "n_ms1_spectra": 0,
            "n_ms2_spectra": 0,
            "fragmentation_type": "unknown",
            "high_accuracy_precursors": "unknown",
        }

        # Just count spectra from file
        if str(mzml_file).endswith(".gz"):
            infile = gzip.open(mzml_file, "rt", encoding="utf-8", errors="ignore")
        else:
            infile = open(mzml_file, "r", encoding="utf-8", errors="ignore")

        try:
            for line in infile:
                if "<spectrum " in line:
                    stats["n_spectra"] += 1
                if 'ms level="1"' in line or "ms level=1" in line:
                    stats["n_ms1_spectra"] += 1
                if 'ms level="2"' in line or "ms level=2" in line:
                    stats["n_ms2_spectra"] += 1
                if "@hcd" in line.lower():
                    stats["fragmentation_type"] = "HR_HCD"
                if "@cid" in line.lower():
                    stats["fragmentation_type"] = "LR_IT_CID"
        finally:
            infile.close()

        return stats

    def _parse_filter_string(self, filter_string: str, stats: Dict) -> None:
        """
        Parse Thermo filter string for fragmentation type.

        Adapted from RunAssessor's parse_filter_string method.
        Handles both MS2 and MS3 fragmentation types separately.

        Args:
            filter_string: Thermo filter string
            stats: Statistics dictionary to update
        """
        if not filter_string:
            return

        # Determine mass accuracy from detector
        mass_accuracy = "??"
        match = re.search(r"^(\S+)", filter_string)
        if match:
            detector = match.group(0)
            if detector in ["FTMS", "ASTMS"]:
                mass_accuracy = "HR"
            elif detector == "ITMS":
                mass_accuracy = "LR"

        # Determine MS level
        ms_level = 0
        if " ms " in filter_string:
            ms_level = 1
        else:
            match = re.search(r" ms(\d) ", filter_string)
            if match:
                ms_level = int(match.group(1))

        # Check fragmentation types
        have_hcd = "@hcd" in filter_string.lower()
        have_cid = "@cid" in filter_string.lower()
        have_etd = "@etd" in filter_string.lower()

        if ms_level == 1:
            stats["n_ms1_spectra"] += 1
            if mass_accuracy == "HR":
                if stats["high_accuracy_precursors"] == "unknown":
                    stats["high_accuracy_precursors"] = "true"
            elif mass_accuracy == "LR":
                if stats["high_accuracy_precursors"] == "unknown":
                    stats["high_accuracy_precursors"] = "false"

        elif ms_level == 2:
            stats["n_ms2_spectra"] += 1
            self._update_fragmentation_stats(stats, "ms2", mass_accuracy, have_hcd, have_cid, have_etd)

        elif ms_level == 3:
            stats["n_ms3_spectra"] += 1
            self._update_fragmentation_stats(stats, "ms3", mass_accuracy, have_hcd, have_cid, have_etd)

    def _update_fragmentation_stats(
        self, stats: Dict, ms_level: str, mass_accuracy: str,
        have_hcd: bool, have_cid: bool, have_etd: bool
    ) -> None:
        """Update fragmentation stats for a given MS level."""
        dissociation_sum = int(have_hcd) + int(have_cid) + int(have_etd)

        if dissociation_sum == 1:
            if have_hcd:
                stats[f"n_{ms_level}_{mass_accuracy}_HCD_spectra"] += 1
            elif have_cid:
                stats[f"n_{ms_level}_{mass_accuracy}_IT_CID_spectra"] += 1
            elif have_etd:
                stats[f"n_{ms_level}_{mass_accuracy}_IT_ETD_spectra"] += 1
        elif dissociation_sum == 2:
            if have_etd and have_hcd:
                stats[f"n_{ms_level}_{mass_accuracy}_EThcD_spectra"] += 1
            elif have_etd and have_cid:
                stats[f"n_{ms_level}_{mass_accuracy}_ETciD_spectra"] += 1

    def _detect_fragmentation_cv(self, spectrum: Dict, stats: Dict) -> None:
        """Detect fragmentation type from CV terms."""
        # Check for dissociation method CV terms
        if "HCD" in str(spectrum) or "beam-type collision-induced dissociation" in str(
            spectrum
        ):
            stats["n_ms2_HR_HCD_spectra"] = stats.get("n_ms2_HR_HCD_spectra", 0) + 1
        elif "CID" in str(spectrum) or "collision-induced dissociation" in str(
            spectrum
        ):
            stats["n_ms2_LR_IT_CID_spectra"] = stats.get("n_ms2_LR_IT_CID_spectra", 0) + 1

    def _determine_fragmentation_type(self, stats: Dict) -> None:
        """Determine dominant fragmentation type for MS2 and MS3 separately."""
        # Determine MS2 fragmentation type
        stats["ms2_fragmentation_type"] = self._get_dominant_fragmentation(stats, "ms2")

        # Determine MS3 fragmentation type
        stats["ms3_fragmentation_type"] = self._get_dominant_fragmentation(stats, "ms3")

        # Set primary fragmentation_type (prefer MS2, fall back to MS3)
        if stats["ms2_fragmentation_type"] != "unknown":
            stats["fragmentation_type"] = stats["ms2_fragmentation_type"]
        elif stats["ms3_fragmentation_type"] != "unknown":
            stats["fragmentation_type"] = stats["ms3_fragmentation_type"]

    def _get_dominant_fragmentation(self, stats: Dict, ms_level: str) -> str:
        """Get the dominant fragmentation type for a given MS level."""
        frag_counts = {}
        for prefix in ["HR", "LR"]:
            for frag in ["HCD", "IT_CID", "IT_ETD", "EThcD", "ETciD"]:
                key = f"n_{ms_level}_{prefix}_{frag}_spectra"
                count = stats.get(key, 0)
                if count > 0:
                    frag_counts[f"{prefix}_{frag}"] = count

        if not frag_counts:
            return "unknown"

        dominant = max(frag_counts, key=frag_counts.get)
        if len(frag_counts) > 1:
            # Check if there's really multiple types or just one dominant
            total = sum(frag_counts.values())
            if frag_counts[dominant] / total < 0.9:
                return "multiple"
        return dominant

    def _estimate_tolerances(self, result: AnalysisResult) -> None:
        """
        Estimate mass tolerances based on instrument category.

        Args:
            result: AnalysisResult to update
        """
        category = result.instrument_category
        frag_type = result.fragmentation_type

        # Default tolerances based on instrument category
        if category in ["pureHCD", "variable"]:
            # High resolution instruments
            result.precursor_tolerance_ppm = 10.0
            if "HR" in frag_type or category == "pureHCD":
                result.fragment_tolerance_ppm = 20.0
                result.confidence["fragment_tolerance"] = 0.7
            else:
                result.fragment_tolerance_da = 0.6
                result.confidence["fragment_tolerance"] = 0.6
            result.confidence["precursor_tolerance"] = 0.7

        elif category == "ion_trap":
            # Low resolution instruments
            result.precursor_tolerance_ppm = 500.0
            result.fragment_tolerance_da = 0.6
            result.confidence["precursor_tolerance"] = 0.6
            result.confidence["fragment_tolerance"] = 0.6

        elif category == "QTOF":
            # Q-TOF instruments
            result.precursor_tolerance_ppm = 20.0
            result.fragment_tolerance_ppm = 40.0
            result.confidence["precursor_tolerance"] = 0.6
            result.confidence["fragment_tolerance"] = 0.6

        else:
            # Unknown - use conservative defaults
            result.precursor_tolerance_ppm = 20.0
            result.fragment_tolerance_ppm = 20.0
            result.confidence["precursor_tolerance"] = 0.4
            result.confidence["fragment_tolerance"] = 0.4

    def analyze_multiple(self, mzml_files: List[Path]) -> List[AnalysisResult]:
        """
        Analyze multiple mzML files.

        Args:
            mzml_files: List of paths to mzML files

        Returns:
            List of analysis results
        """
        results = []
        for mzml_file in mzml_files:
            result = self.analyze(mzml_file)
            results.append(result)
        return results

    def aggregate_results(self, results: List[AnalysisResult]) -> Dict[str, Any]:
        """
        Aggregate results from multiple files into consensus values.

        Args:
            results: List of analysis results

        Returns:
            Aggregated/consensus parameters
        """
        if not results:
            return {}

        aggregated: Dict[str, Any] = {
            "instrument_model": {},
            "fragmentation_type": {},
            "ms2_fragmentation_type": {},
            "ms3_fragmentation_type": {},
            "precursor_tolerance_ppm": [],
            "fragment_tolerance_ppm": [],
            "fragment_tolerance_da": [],
            "high_accuracy_precursors": [],
            "acquisition_type": {},
            "ms_min_charge": [],
            "ms_max_charge": [],
            "collision_energy": {},
            # pyopenms-derived fields
            "mass_analyzer_type": {},
            "ionization_method": {},
            "polarity": {},
            "ms1_scan_range": [],
            "isolation_window_da": [],
            "instrument_serial_number": [],
        }

        for r in results:
            # Count instrument models
            inst_name = r.instrument_model.get("name", "unknown")
            aggregated["instrument_model"][inst_name] = (
                aggregated["instrument_model"].get(inst_name, 0) + 1
            )

            # Count fragmentation types (overall, MS2, MS3)
            frag = r.fragmentation_type
            aggregated["fragmentation_type"][frag] = (
                aggregated["fragmentation_type"].get(frag, 0) + 1
            )

            ms2_frag = r.ms2_fragmentation_type
            if ms2_frag != "unknown":
                aggregated["ms2_fragmentation_type"][ms2_frag] = (
                    aggregated["ms2_fragmentation_type"].get(ms2_frag, 0) + 1
                )

            ms3_frag = r.ms3_fragmentation_type
            if ms3_frag != "unknown":
                aggregated["ms3_fragmentation_type"][ms3_frag] = (
                    aggregated["ms3_fragmentation_type"].get(ms3_frag, 0) + 1
                )

            # Collect tolerances
            if r.precursor_tolerance_ppm is not None:
                aggregated["precursor_tolerance_ppm"].append(r.precursor_tolerance_ppm)
            if r.fragment_tolerance_ppm is not None:
                aggregated["fragment_tolerance_ppm"].append(r.fragment_tolerance_ppm)
            if r.fragment_tolerance_da is not None:
                aggregated["fragment_tolerance_da"].append(r.fragment_tolerance_da)

            aggregated["high_accuracy_precursors"].append(r.high_accuracy_precursors)

            acq = r.acquisition_type
            aggregated["acquisition_type"][acq] = (
                aggregated["acquisition_type"].get(acq, 0) + 1
            )

            # Collect charge states
            if r.ms_min_charge is not None:
                aggregated["ms_min_charge"].append(r.ms_min_charge)
            if r.ms_max_charge is not None:
                aggregated["ms_max_charge"].append(r.ms_max_charge)

            # Collect collision energy
            if r.collision_energy:
                aggregated["collision_energy"][r.collision_energy] = (
                    aggregated["collision_energy"].get(r.collision_energy, 0) + 1
                )

            # Collect pyopenms-derived fields
            if r.mass_analyzer_type != "unknown":
                aggregated["mass_analyzer_type"][r.mass_analyzer_type] = (
                    aggregated["mass_analyzer_type"].get(r.mass_analyzer_type, 0) + 1
                )
            if r.ionization_method != "unknown":
                aggregated["ionization_method"][r.ionization_method] = (
                    aggregated["ionization_method"].get(r.ionization_method, 0) + 1
                )
            if r.polarity != "unknown":
                aggregated["polarity"][r.polarity] = (
                    aggregated["polarity"].get(r.polarity, 0) + 1
                )
            if r.ms1_scan_range is not None:
                aggregated["ms1_scan_range"].append(r.ms1_scan_range)
            if r.isolation_window_da is not None:
                aggregated["isolation_window_da"].append(r.isolation_window_da)
            if r.instrument_serial_number:
                aggregated["instrument_serial_number"].append(r.instrument_serial_number)

        # Compute consensus values
        consensus = {}

        # Most common instrument
        if aggregated["instrument_model"]:
            most_common_name = max(
                aggregated["instrument_model"], key=aggregated["instrument_model"].get
            )
            # Get full instrument data from first matching result
            instrument_dict = None
            for r in results:
                if r.instrument_model.get("name") == most_common_name:
                    instrument_dict = r.instrument_model.copy()
                    break
            # If no matching dict found, create one from the name
            if instrument_dict:
                consensus["instrument_model"] = instrument_dict
            else:
                consensus["instrument_model"] = {"name": most_common_name}
            consensus["instrument_model_confidence"] = (
                aggregated["instrument_model"].get(most_common_name, 0) / len(results)
            )

        # Most common fragmentation (overall)
        if aggregated["fragmentation_type"]:
            consensus["fragmentation_type"] = max(
                aggregated["fragmentation_type"],
                key=aggregated["fragmentation_type"].get,
            )

        # Most common MS2 fragmentation
        if aggregated["ms2_fragmentation_type"]:
            consensus["ms2_fragmentation_type"] = max(
                aggregated["ms2_fragmentation_type"],
                key=aggregated["ms2_fragmentation_type"].get,
            )

        # Most common MS3 fragmentation
        if aggregated["ms3_fragmentation_type"]:
            consensus["ms3_fragmentation_type"] = max(
                aggregated["ms3_fragmentation_type"],
                key=aggregated["ms3_fragmentation_type"].get,
            )

        # Median tolerances (more robust than mean)
        if aggregated["precursor_tolerance_ppm"]:
            consensus["precursor_tolerance_ppm"] = float(
                np.median(aggregated["precursor_tolerance_ppm"])
            )
        if aggregated["fragment_tolerance_ppm"]:
            consensus["fragment_tolerance_ppm"] = float(
                np.median(aggregated["fragment_tolerance_ppm"])
            )
        if aggregated["fragment_tolerance_da"]:
            consensus["fragment_tolerance_da"] = float(
                np.median(aggregated["fragment_tolerance_da"])
            )

        # High accuracy consensus
        consensus["high_accuracy_precursors"] = (
            sum(aggregated["high_accuracy_precursors"]) > len(results) / 2
        )

        # Acquisition type
        if aggregated["acquisition_type"]:
            consensus["acquisition_type"] = max(
                aggregated["acquisition_type"],
                key=aggregated["acquisition_type"].get,
            )

        # Charge state range (use min of mins and max of maxes across files)
        if aggregated["ms_min_charge"]:
            consensus["ms_min_charge"] = min(aggregated["ms_min_charge"])
        if aggregated["ms_max_charge"]:
            consensus["ms_max_charge"] = max(aggregated["ms_max_charge"])

        # Collision energy (most common across files)
        if aggregated["collision_energy"]:
            consensus["collision_energy"] = max(
                aggregated["collision_energy"],
                key=aggregated["collision_energy"].get,
            )

        # pyopenms-derived consensus
        if aggregated["mass_analyzer_type"]:
            raw_ma = max(
                aggregated["mass_analyzer_type"],
                key=aggregated["mass_analyzer_type"].get,
            )
            ont = MASS_ANALYZER_ONTOLOGY.get(raw_ma)
            if ont:
                consensus["mass_analyzer_type"] = {
                    "name": ont["name"],
                    "accession": ont["accession"],
                    "raw": raw_ma,
                }
            else:
                consensus["mass_analyzer_type"] = {"name": raw_ma, "accession": "", "raw": raw_ma}

        if aggregated["ionization_method"]:
            raw_ion = max(
                aggregated["ionization_method"],
                key=aggregated["ionization_method"].get,
            )
            ont = IONIZATION_ONTOLOGY.get(raw_ion)
            if ont:
                consensus["ionization_method"] = {
                    "name": ont["name"],
                    "accession": ont["accession"],
                    "raw": raw_ion,
                }
            else:
                consensus["ionization_method"] = {"name": raw_ion, "accession": "", "raw": raw_ion}

        if aggregated["polarity"]:
            pol_keys = set(aggregated["polarity"].keys())
            if pol_keys == {"positive"}:
                consensus["polarity"] = "positive"
            elif pol_keys == {"negative"}:
                consensus["polarity"] = "negative"
            elif len(pol_keys) > 1 or "mixed" in pol_keys:
                consensus["polarity"] = "mixed"
            else:
                consensus["polarity"] = next(iter(pol_keys))

        if aggregated["ms1_scan_range"]:
            # Use first (they should be consistent across files)
            consensus["ms1_scan_range"] = aggregated["ms1_scan_range"][0]
        if aggregated["isolation_window_da"]:
            consensus["isolation_window_da"] = float(
                np.median(aggregated["isolation_window_da"])
            )
        if aggregated["instrument_serial_number"]:
            # Use first non-empty
            consensus["instrument_serial_number"] = aggregated["instrument_serial_number"][0]

        # Annotate acquisition type with ontology
        acq_raw = consensus.get("acquisition_type")
        if acq_raw and acq_raw in ACQUISITION_ONTOLOGY:
            ont = ACQUISITION_ONTOLOGY[acq_raw]
            consensus["acquisition_type_ont"] = {
                "name": ont["name"],
                "accession": ont["accession"],
                "raw": acq_raw,
            }

        return consensus
