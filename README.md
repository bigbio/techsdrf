# TechSDRF

A tool to validate and refine SDRF (Sample and Data Relationship Format) files using empirical mass spectrometry data analysis.

## Overview

TechSDRF downloads raw MS files from PRIDE or MassIVE, analyzes them to detect instrument parameters and post-translational modifications, and proposes refinements to SDRF metadata files. It helps ensure that SDRF annotations accurately reflect the actual experimental parameters.

## Features

- **Automatic parameter detection** from raw MS data:
  - Instrument model and serial number
  - Fragmentation method (HCD, CID, ETD, EThcD) with MS2/MS3 support
  - Precursor and fragment mass tolerances (empirical Gaussian estimation)
  - Acquisition type (DDA/DIA)
  - Charge state range, collision energy, mass analyzer, ionization method
  - Scan window range, isolation window width, polarity

- **PTM detection** from spectra without database search:
  - Isobaric labels (TMT6/11/16/18plex, iTRAQ4/8plex) via reporter ion detection
  - Enrichment signatures (phospho, glyco) via diagnostic ions
  - Common variable modifications via precursor mass pairing with statistical scoring
  - See [Algorithms](#algorithms) for details

- **Multi-repository support**:
  - PRIDE (PXD accessions) via PRIDE API
  - MassIVE (MSV accessions) via FTP

- **Multi-vendor support**:
  - Thermo RAW files (via ThermoRawFileParser)
  - Bruker .d files (via msconvert)
  - AB SCIEX .wiff files (via msconvert)

- **Comprehensive reporting**:
  - Plain text report with recommendations and confidence levels
  - PTM detection results with statistical scores
  - Refined SDRF output

## Installation

### Using Conda (Recommended)

```bash
conda env create -f environment.yml
conda activate techsdrf
pip install -e .
```

### Using pip

```bash
pip install techsdrf
```

**Note**: ThermoRawFileParser must be installed separately for Thermo RAW file support:
```bash
conda install -c bioconda thermorawfileparser
```

## Usage

### Basic Usage

```bash
# Refine an SDRF file using data from PRIDE
techsdrf refine -p PXD012345 -s input.sdrf.tsv

# From MassIVE
techsdrf refine -p MSV000085836 -s input.sdrf.tsv

# With more files for higher confidence
techsdrf refine -p PXD012345 -s input.sdrf.tsv -n 5

# Analyze all files
techsdrf refine -p PXD012345 -s input.sdrf.tsv --all-files

# Skip PTM detection
techsdrf refine -p PXD012345 -s input.sdrf.tsv --skip-ptm

# Verbose output
techsdrf refine -p PXD012345 -s input.sdrf.tsv -vv
```

### Options

```
techsdrf refine [OPTIONS]

Options:
  -p, --accession TEXT              Dataset accession (PXD or MSV) [required]
  -s, --sdrf-file PATH              Path to existing SDRF file [required]
  -n, --num-files INTEGER           Number of raw files to analyze (default: 3)
  -a, --all-files                   Analyze all raw files
  -o, --output-sdrf PATH            Output path for refined SDRF
  -r, --report-file PATH            Output path for report
  -w, --work-dir PATH               Working directory for downloads
  -t, --file-types [raw|d|wiff|wiff2]  Filter by file type
  -u, --tolerance-unit [ppm|Da]     Preferred unit for tolerances
  --keep-files                      Keep downloaded/converted files
  --skip-ptm                        Skip PTM detection from spectra
  -v, --verbose                     Increase verbosity (-v, -vv)
```

### Other Commands

```bash
techsdrf validate -s input.sdrf.tsv   # Validate an SDRF file
techsdrf info -s input.sdrf.tsv       # Show SDRF parameters
```

## Algorithms

### Instrument Parameter Detection

TechSDRF extracts instrument parameters from mzML files using header parsing for static metadata and spectrum-level scanning for runtime parameters.

**Header parsing (instrument model, mass analyzer, ionization):**
1. Load the mzML file via **pyopenms** (`MSExperiment` / `MzMLFile`).
2. Match the instrument name/model string against a curated lookup table of ~40 Thermo, Bruker, Waters, and AB SCIEX instruments, each mapped to a PSI-MS ontology accession (e.g., `MS:1002416` for Orbitrap Fusion) and an instrument category (`pureHCD`, `variable`, `ion_trap`, `QTOF`).
3. Extract mass analyzer type, ionization method, and instrument serial number from the pyopenms `Instrument` object.

**Spectrum-level scanning (fragmentation, charge, acquisition, polarity):**
1. Iterate up to 10,000 spectra. For Thermo data, parse the **filter string** (e.g., `FTMS + p NSI Full ms2 400.00@hcd30.00`) to extract:
   - Detector type -> mass accuracy (`FTMS`/`ASTMS` = high resolution, `ITMS` = low resolution)
   - Fragmentation method -> `HCD`, `CID`, `ETD`, `EThcD`, `ETciD` (supports both MS2 and MS3)
   - MS level (ms1, ms2, ms3)
2. For non-Thermo data, use CV terms from the mzML spectrum metadata.
3. **Fragmentation consensus:** Count spectra per fragmentation type per MS level. The dominant type (>90% of spectra) is reported; otherwise `multiple`.
4. **Acquisition type:** Determined from median isolation window width: <=3 Da -> DDA, >=15 Da -> DIA.
5. **Charge state range:** Collected from precursor `charge state` CV terms across all MS2 spectra.
6. **Collision energy:** Extracted from `MS:1000045` (collision energy) and `MS:1000138` (percent collision energy) activation parameters. Stepped CE is detected when 2+ distinct values each appear in >5% of spectra.
7. **Polarity, scan window, isolation window:** Aggregated from pyopenms spectrum-level instrument settings.

### Tolerance Estimation

Precursor and fragment mass tolerances are estimated empirically using two independent methods, with a consensus algorithm to select the best estimate.

**Method 1: Gaussian fitting (RunAssessor)**
- Uses the [RunAssessor](https://github.com/mpc-bioinformatics/RunAssessor) library to fit Gaussian distributions to mass error profiles.
- RunAssessor reads spectra from the mzML file, builds composite spectra from low-end and neutral-loss regions, identifies regions of interest (ROIs), and fits the mass error distribution to estimate the standard deviation.
- Reports precursor tolerance in ppm and fragment tolerance in ppm or Da (depending on the mass analyzer type).

**Method 2: Crux param-medic**
- Uses the [Crux](http://crux.ms/) `param-medic` tool, which estimates mass measurement accuracy from paired precursor/fragment peak distributions.
- Runs as an external subprocess with a 5-minute timeout.
- Reports precursor error in ppm and fragment bin size in Th (Da).
- Note: Requires a compatibility patch for mzML files from ThermoRawFileParser (PSI term `1003145` -> `1000615`), applied automatically to a temporary copy.

**Consensus algorithm:**
1. Collect estimates from both methods (when available).
2. Check agreement: if all estimates are within 20% of the median, compute a **confidence-weighted average**.
3. If a majority agrees (>=2 methods within 20%), use the **median** (robust to outliers).
4. If methods disagree significantly, select the **highest-confidence** estimate, with Gaussian preferred over Crux when confidence is tied.
5. The consensus reports precursor and fragment tolerances in both ppm and Da, using the native unit from the winning method.

**Instrument-based heuristic fallback:**
When empirical estimation is unavailable, default tolerances are assigned by instrument category:

| Category | Precursor | Fragment |
|----------|-----------|----------|
| pureHCD (Exactive, Q Exactive) | 10 ppm | 20 ppm |
| variable (Fusion, Lumos, Eclipse) - HR frag | 10 ppm | 20 ppm |
| variable - LR frag (IT-CID) | 10 ppm | 0.6 Da |
| ion_trap (LTQ, Velos) | 500 ppm | 0.6 Da |
| QTOF (timsTOF, TripleTOF) | 20 ppm | 40 ppm |

### PTM Detection

TechSDRF detects post-translational modifications directly from MS spectra using three independent detection tiers. Results are reported but not automatically written to SDRF modification columns.

#### Tier 1: Reporter Ion Detection

Detects isobaric labeling reagents (TMT, iTRAQ) by scanning the low-m/z region (100–140 m/z) of MS2/MS3 spectra for characteristic reporter ion clusters.

**Method:**
1. For each MS2/MS3 spectrum, extract peaks in the 100–140 m/z window.
2. Match peaks against known reporter ion m/z values (±0.01 Da tolerance).
3. Count how many channels are present per spectrum. A spectrum is considered labeled if the number of matched channels meets the minimum threshold (e.g., ≥4 of 6 for TMT6plex).
4. Count the fraction of spectra with reporter ions to determine overall confidence.

**Plex-level disambiguation:**
TMT plex variants (6/11/16/18) share reporter channels — TMT6 channels are a subset of TMT11, which is a subset of TMT16, etc. To distinguish plex levels:
- For each higher plex, identify channels *unique* to that plex (not present in the lower set).
- Verify unique channels with tight tolerance (0.003 Da) and ensure each matched peak is closer to the unique channel than to any base-set channel.
- Require that at least half of the unique channels are confirmed in ≥50% of labeled spectra.
- Iterate from most specific (TMT18) to least specific (TMT6), reporting the first plex that passes validation.

**Confidence scoring:**
| Fraction of spectra | Confidence |
|---------------------|------------|
| ≥50%                | HIGH (0.95)|
| ≥20%                | MEDIUM (0.85)|
| ≥10%                | MEDIUM (0.70)|
| ≥5%                 | LOW (0.50) |

#### Tier 2: Diagnostic Ion Detection

Identifies PTMs that produce characteristic fragment signatures: neutral losses (phosphorylation) and oxonium ions (glycosylation).

**Neutral losses (phosphorylation):**
- For each MS2 spectrum with a known precursor, compute the expected neutral loss peak position: `target_mz = precursor_mz − neutral_loss_Da / charge`.
- Check if a peak exists at this position (±0.02 Da) with relative intensity ≥5% of the base peak.
- Two neutral losses are checked: H₃PO₄ (97.977 Da) and HPO₃ (79.966 Da).

**Oxonium ions (glycosylation):**
- Check for the presence of HexNAc oxonium ion at 204.087 m/z and Hex(1)HexNAc(1) at 366.140 m/z.
- Must have relative intensity ≥5% (HexNAc) or ≥2% (HexHexNAc) of the base peak.

**Confidence scoring:**
| Fraction of spectra | Confidence |
|---------------------|------------|
| ≥20%                | HIGH (0.90)|
| ≥10%                | MEDIUM (0.70)|
| ≥5%                 | LOW (0.50) |
| ≥1%                 | LOW (0.30) |

#### Tier 3: Open Mass-Shift Detection

Detects common variable modifications by looking for precursor mass pairs separated by known PTM mass shifts, with statistical scoring to distinguish real modifications from random coincidences.

**Method:**
1. **Collect and filter precursor neutral masses** from MS2 spectra:
   `neutral_mass = (precursor_mz − 1.00728) × charge`

   Quality gates applied before mass pairing:
   - Minimum 15 fragment peaks (removes noise triggers)
   - Charge state ≥ 2 (skips unknown and singly-charged non-peptide precursors)
   - Neutral mass in 400–6000 Da (peptide range)

2. **Deduplicate to unique species:** Masses are binned to 0.01 Da resolution and deduplicated. This converts the metric from "number of MS2 scans" to "number of distinct precursor species", reducing inflation from repeated fragmentation of the same peptide and improving the accuracy of the null model which assumes uniform mass density.

3. **Count mass pairs** for each known PTM shift: for every unique precursor mass *M* in the sorted array, use binary search to check whether *M* + Δ (the PTM mass shift) exists within ±0.02 Da. This identifies unmodified/modified peptide pairs.

4. **Null model — expected random matches:** Under the assumption that unique precursor masses are uniformly distributed:
   ```
   density = N / mass_range          (unique masses per Da)
   E[matches] = N × density × 2×tol  (expected random pairs)
   ```

5. **Enrichment scoring:**
   ```
   enrichment = observed_count / expected_count
   ```
   An enrichment of 4× means the modification is observed 4 times more often than expected by chance.

6. **Statistical significance** via Poisson survival function:
   ```
   p-value = P(X ≥ observed | λ = expected)
   probability = 1 − p-value
   ```
   Uses `scipy.stats.poisson.sf` when available, with fallback to normal approximation (large λ) or direct log-space summation (small λ).

7. **Filtering and ranking:**
   - Enrichment > 1.5×
   - p-value < 0.01
   - Minimum 5 unique pairs
   - Ranked by enrichment, top 10 reported

**Known mass shifts searched:**

| Modification | Δ mass (Da) | UNIMOD |
|-------------|-------------|--------|
| Oxidation | +15.995 | UNIMOD:35 |
| Deamidation | +0.984 | UNIMOD:7 |
| Phospho | +79.966 | UNIMOD:21 |
| Acetyl | +42.011 | UNIMOD:1 |
| Methyl | +14.016 | UNIMOD:34 |
| Dimethyl | +28.031 | UNIMOD:36 |
| Carbamidomethyl | +57.022 | UNIMOD:4 |
| Carbamyl | +43.006 | UNIMOD:5 |
| Propionamide | +71.037 | UNIMOD:24 |
| Pyro-glu from Q | −17.027 | UNIMOD:28 |
| Pyro-glu from E | −18.011 | UNIMOD:27 |
| Formyl | +27.995 | UNIMOD:122 |
| Dioxidation | +31.990 | UNIMOD:425 |

**Confidence scoring (from enrichment):**
| Enrichment | Confidence |
|-----------|------------|
| ≥5.0× | HIGH (0.95) |
| ≥3.0× | MEDIUM (0.85) |
| ≥2.0× | MEDIUM (0.70) |
| ≥1.5× | LOW (0.50) |

#### Fixed vs Variable Classification

Each detected PTM is classified as **fixed** or **variable** based on its tier and identity:

| Tier | Rule | Examples |
|------|------|----------|
| Tier 1 (reporter ions) | Always **fixed** | TMT, iTRAQ labels are covalently attached to all peptides |
| Tier 2 (diagnostic ions) | Always **variable** | Phosphorylation, glycosylation occur on a subset of peptides |
| Tier 3 (mass shifts) | Lookup by UNIMOD accession | Carbamidomethyl (UNIMOD:4) and Propionamide (UNIMOD:24) are **fixed** (alkylation reagents applied to all cysteines); all others are **variable** |

#### Multi-File Aggregation

When multiple files are analyzed, per-file PTM results are aggregated:
- Evidence counts and spectra scanned are summed across files.
- The maximum confidence, enrichment, and probability across files are reported.
- Reporter ion type is determined by majority vote.

## Example Report

```
PTM DETECTION
----------------------------------------
  Files analyzed for PTMs: 3

  Reporter Ions:
    TMT11plex (UNIMOD:737): detected in 52% of spectra (3/3 files) [HIGH] [fixed]

  Diagnostic Ions:
    Phospho (UNIMOD:21): 8.3% of spectra (415 observations) [MEDIUM] [variable]

  Mass Shifts:
    Oxidation (+15.995 Da) (UNIMOD:35): 2840 obs (3/3 files) [MEDIUM] [variable] (enrichment=4.2x, prob=1.00)
    Deamidation (+0.984 Da) (UNIMOD:7): 2650 obs (3/3 files) [MEDIUM] [variable] (enrichment=3.9x, prob=1.00)
    Methyl (+14.016 Da) (UNIMOD:34): 2410 obs (3/3 files) [MEDIUM] [variable] (enrichment=3.6x, prob=1.00)
```

## Possible Improvements

- **SAGE integration**: Empirical tolerance estimation via [Sage](https://github.com/lazear/sage) database search optimization could be added to refine precursor and fragment mass tolerances by grid search over a FASTA database. This would complement the existing Gaussian (RunAssessor) and Crux param-medic methods.
- **PTM auto-apply**: Optional `--apply-ptm` flag to write detected PTMs to SDRF `comment[modification parameters]` columns.

## Dependencies

- Python >= 3.9
- numpy, scipy, pandas
- pyopenms (for mzML parsing and spectrum access)
- sdrf-pipelines (for SDRF validation)
- pridepy (for PRIDE API access)
- runassessor (for Gaussian tolerance estimation)
- ThermoRawFileParser (for Thermo RAW conversion)
- ProteoWizard msconvert (for Bruker/SCIEX conversion)

## License

Apache 2.0
