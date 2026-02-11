"""
Configuration constants and defaults for SDRF Refiner.

Ontology catalog
~~~~~~~~~~~~~~~~
All CV term mappings in this module follow SDRF format conventions:

* Ontology-valued columns use ``AC=<accession>;NT=<name>``  (or the reverse
  ``NT=…;AC=…``).  Both orderings are accepted by sdrf-pipelines.
* Tolerance columns use ``<value> <unit>``  (e.g. ``10 ppm``, ``0.02 Da``).
* Plain-text columns (charge, collision energy) use raw values.

The canonical ontology sources are:

* **PSI-MS** – mass spectrometry ontology
  https://www.ebi.ac.uk/ols4/ontologies/ms
* **PRIDE** – PRIDE controlled vocabulary
  https://www.ebi.ac.uk/ols4/ontologies/pride
* **UNIMOD** – protein modifications
  https://www.unimod.org
"""

from typing import Any, Dict

# =============================================================================
# General settings
# =============================================================================

# Default number of files to download from PRIDE
DEFAULT_NUM_FILES = 3

# Supported raw file extensions by vendor
SUPPORTED_FILE_TYPES = ["raw", "d", "wiff", "wiff2"]

# Vendor to extension mapping
VENDOR_EXTENSIONS = {
    "thermo": [".raw"],
    "bruker": [".d"],
    "sciex": [".wiff", ".wiff2"],
    "waters": [".raw"],  # Waters also uses .raw
}

# PRIDE API base URL
PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v2"

# Tolerance comparison threshold (20% difference is acceptable)
TOLERANCE_MATCH_THRESHOLD = 0.2

# Soft divergence limit: when a single estimation method (or disagreeing
# methods) produces a value that differs from the SDRF by more than this
# factor, the change is flagged for review.  Bypassed when multiple methods
# agree (strong consensus).
TOLERANCE_DIVERGENCE_SOFT = 5.0

# Hard divergence limit: when the empirical value differs from the SDRF by
# more than this factor, the change is ALWAYS flagged for review — even if
# multiple methods agree.  A 10x+ difference likely means the SDRF captures
# a search-engine parameter while the empirical value reflects instrument
# accuracy; these are fundamentally different things.
TOLERANCE_DIVERGENCE_HARD = 10.0

# =============================================================================
# SDRF column names
# =============================================================================
# Maps internal parameter keys → SDRF column headers as defined in the
# sdrf-templates YAML schemas (dda-acquisition, dia-acquisition, ms-proteomics).

SDRF_COLUMNS: Dict[str, str] = {
    "instrument": "comment[instrument]",
    "dissociation": "comment[dissociation method]",
    "ms2_dissociation": "comment[ms2 dissociation method]",
    "ms3_dissociation": "comment[ms3 dissociation method]",
    "precursor_tolerance": "comment[precursor mass tolerance]",
    "fragment_tolerance": "comment[fragment mass tolerance]",
    "data_file": "comment[data file]",
    "acquisition_method": "comment[proteomics data acquisition method]",
    "ms_min_charge": "comment[MS min charge]",
    "ms_max_charge": "comment[MS max charge]",
    "collision_energy": "comment[collision energy]",
    # Columns populated from pyopenms header analysis
    "ms2_mass_analyzer": "comment[ms2 mass analyzer]",
    # DIA-specific (scan window + isolation width)
    "scan_window_lower": "comment[scan window lower limit]",
    "scan_window_upper": "comment[scan window upper limit]",
    "isolation_window_width": "comment[isolation window width]",
}

# =============================================================================
# Dissociation / fragmentation ontology   (PSI-MS children of MS:1000044)
# =============================================================================
# Maps the internal fragmentation-type keys (e.g. "HR_HCD") that come out of
# MSAnalyzer to the canonical PSI-MS term name and accession.

FRAGMENTATION_ONTOLOGY: Dict[str, Dict[str, str]] = {
    "HR_HCD":    {"name": "HCD",   "accession": "MS:1000422"},
    "LR_HCD":    {"name": "HCD",   "accession": "MS:1000422"},
    "HR_IT_CID": {"name": "CID",   "accession": "MS:1000133"},
    "LR_IT_CID": {"name": "CID",   "accession": "MS:1000133"},
    "HR_IT_ETD": {"name": "ETD",   "accession": "MS:1000598"},
    "LR_IT_ETD": {"name": "ETD",   "accession": "MS:1000598"},
    "HR_EThcD":  {"name": "EThcD", "accession": "MS:1002631"},
    "LR_EThcD":  {"name": "EThcD", "accession": "MS:1002631"},
    "HR_ETciD":  {"name": "ETciD", "accession": "MS:1002632"},
    "LR_ETciD":  {"name": "ETciD", "accession": "MS:1002632"},
    "HR_QTOF":   {"name": "CID",   "accession": "MS:1000133"},
}

# =============================================================================
# Mass analyzer ontology   (PSI-MS children of MS:1000443 "mass analyzer type")
# =============================================================================
# Maps the pyopenms MassAnalyzer.AnalyzerType enum name (as produced by
# MSAnalyzer._MASS_ANALYZER_NAMES) to the canonical PSI-MS term.

MASS_ANALYZER_ONTOLOGY: Dict[str, Dict[str, str]] = {
    "Orbitrap":              {"name": "orbitrap",                               "accession": "MS:1000484"},
    "IT":                    {"name": "ion trap",                               "accession": "MS:1000264"},
    "LIT":                   {"name": "linear ion trap",                        "accession": "MS:1000291"},
    "PaulIonTrap":           {"name": "ion trap",                               "accession": "MS:1000264"},
    "RadialEjectionLIT":     {"name": "radial ejection linear ion trap",        "accession": "MS:1000083"},
    "AxialEjectionLIT":      {"name": "axial ejection linear ion trap",         "accession": "MS:1000078"},
    "TOF":                   {"name": "time-of-flight",                         "accession": "MS:1000084"},
    "Quadrupole":            {"name": "quadrupole",                             "accession": "MS:1000081"},
    "FourierTransform":      {"name": "fourier transform ion cyclotron resonance mass spectrometer",
                                                                                "accession": "MS:1000079"},
    "Sector":                {"name": "magnetic sector",                        "accession": "MS:1000080"},
    "Cyclotron":             {"name": "fourier transform ion cyclotron resonance mass spectrometer",
                                                                                "accession": "MS:1000079"},
    "IonStorage":            {"name": "stored waveform inverse fourier transform", "accession": "MS:1000284"},
    "ESA":                   {"name": "electrostatic energy analyzer",          "accession": "MS:1000254"},
}

# =============================================================================
# Ionization method ontology   (PSI-MS children of MS:1000008 "ionization type")
# =============================================================================
# Maps pyopenms IonSource.IonizationMethod enum name → PSI-MS term.

IONIZATION_ONTOLOGY: Dict[str, Dict[str, str]] = {
    "ESI":    {"name": "electrospray ionization",                "accession": "MS:1000073"},
    "NESI":   {"name": "nanoelectrospray",                       "accession": "MS:1000398"},
    "MESI":   {"name": "microelectrospray",                      "accession": "MS:1000397"},
    "MALDI":  {"name": "matrix-assisted laser desorption ionization", "accession": "MS:1000075"},
    "AP_MALDI": {"name": "atmospheric pressure matrix-assisted laser desorption ionization",
                                                                  "accession": "MS:1000239"},
    "APCI":   {"name": "atmospheric pressure chemical ionization", "accession": "MS:1000070"},
    "APPI":   {"name": "atmospheric pressure photoionization",    "accession": "MS:1000382"},
    "EI":     {"name": "electron ionization",                     "accession": "MS:1000389"},
    "CI":     {"name": "chemical ionization",                     "accession": "MS:1000071"},
    "FAB":    {"name": "fast atom bombardment ionization",        "accession": "MS:1000074"},
    "FD":     {"name": "field desorption",                        "accession": "MS:1000257"},
    "FI":     {"name": "field ionization",                        "accession": "MS:1000258"},
    "LD":     {"name": "laser desorption ionization",             "accession": "MS:1000393"},
    "PD":     {"name": "plasma desorption ionization",            "accession": "MS:1000400"},
    "SI":     {"name": "surface ionization",                      "accession": "MS:1000406"},
    "TSP":    {"name": "thermospray ionization",                  "accession": "MS:1000407"},
    "ICP":    {"name": "inductively coupled plasma",              "accession": "MS:1000392"},
    "SELDI":  {"name": "surface enhanced laser desorption ionization", "accession": "MS:1000405"},
    "SEND":   {"name": "surface enhanced neat desorption",        "accession": "MS:1000278"},
    "SALDI":  {"name": "surface-assisted laser desorption ionization", "accession": "MS:1000404"},
}

# =============================================================================
# Acquisition method ontology   (PRIDE CV)
# =============================================================================

ACQUISITION_ONTOLOGY: Dict[str, Dict[str, str]] = {
    "DDA": {"name": "Data-Dependent Acquisition", "accession": "PRIDE:0000427"},
    "DIA": {"name": "Data-Independent Acquisition", "accession": "PRIDE:0000451"},
    "SRM": {"name": "Selected reaction monitoring", "accession": "PRIDE:0000553"},
    "PRM": {"name": "Parallel reaction monitoring", "accession": "PRIDE:0000554"},
}

# =============================================================================
# Polarity   (PSI-MS children of MS:1000465 "scan polarity")
# =============================================================================

POLARITY_ONTOLOGY: Dict[str, Dict[str, str]] = {
    "positive": {"name": "positive scan", "accession": "MS:1000130"},
    "negative": {"name": "negative scan", "accession": "MS:1000129"},
}

# =============================================================================
# Default tolerances when detection fails
# =============================================================================

DEFAULT_TOLERANCES: Dict[str, Any] = {
    "HR": {
        "ppm": {"default": 20, "minimum": 5},
        "mz": {"default": 0.01, "minimum": 0.005},
    },
    "LR": {
        "ppm": {"default": 3000, "minimum": 500},
        "mz": {"default": 0.6, "minimum": 0.1},
    },
}


# =============================================================================
# Helpers
# =============================================================================

def format_ontology_value(name: str, accession: str) -> str:
    """Format an ontology value for SDRF output: ``AC=<acc>;NT=<name>``."""
    if name and accession:
        return f"AC={accession};NT={name}"
    return name or ""
