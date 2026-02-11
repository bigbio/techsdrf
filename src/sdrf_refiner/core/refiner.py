"""
SDRFRefiner - Main orchestration class for SDRF refinement workflow.
"""

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional

from sdrf_refiner.analyzer.ms_analyzer import AnalysisResult, MSAnalyzer
from sdrf_refiner.analyzer.ptm_detector import PTMDetectionResult, PTMDetector
from sdrf_refiner.analyzer.tolerance_estimator import ToleranceEstimator
from sdrf_refiner.converter.converter import RawFileConverter
from sdrf_refiner.downloader_factory import create_downloader
from sdrf_refiner.report.generator import ReportGenerator
from sdrf_refiner.sdrf.comparer import ParameterComparer, ParameterComparison
from sdrf_refiner.sdrf.reader import SDRFReader
from sdrf_refiner.sdrf.writer import SDRFWriter
from sdrf_refiner.utils.temp_manager import TempManager

logger = logging.getLogger(__name__)


@dataclass
class RefinementResult:
    """Result of the refinement process."""

    success: bool
    output_sdrf: Optional[Path] = None
    report_file: Optional[Path] = None
    error_message: Optional[str] = None
    comparisons: Optional[List[ParameterComparison]] = None
    detected_params: Optional[Dict[str, Any]] = None
    num_refinements: int = 0


class SDRFRefiner:
    """
    Main class that orchestrates the SDRF refinement workflow.

    Workflow:
    1. Validate input SDRF
    2. Download raw files from PRIDE
    3. Convert raw files to mzML
    4. Analyze mzML files
    5. Compare detected vs SDRF parameters
    6. Generate refinement report
    7. Output refined SDRF
    """

    def __init__(
        self,
        accession: str,
        sdrf_file: Path,
        num_files: Optional[int] = 3,
        output_sdrf: Optional[Path] = None,
        report_file: Optional[Path] = None,
        work_dir: Optional[Path] = None,
        keep_files: bool = False,
        file_types: Optional[List[str]] = None,
        tolerance_unit: Optional[str] = None,
        skip_ptm: bool = False,
    ):
        """
        Initialize SDRF refiner.

        Args:
            accession: Dataset accession (e.g., PXD012345 for PRIDE, MSV000085836 for MassIVE)
            sdrf_file: Path to input SDRF file
            num_files: Number of raw files to download and analyze (None = all files)
            output_sdrf: Output path for refined SDRF (default: <input>_refined.sdrf.tsv)
            report_file: Output path for report (default: <input>_report.txt)
            work_dir: Working directory for downloads and processing (default: cwd)
            keep_files: If True, keep temp files after processing
            file_types: Optional filter for raw file types (e.g., ['raw', 'd'])
            tolerance_unit: Preferred output unit for tolerances ('ppm' or 'Da'); None = auto
            skip_ptm: If True, skip PTM detection from spectra
        """
        self.accession = accession.upper()
        self.pride_accession = self.accession  # backward-compat alias
        self.sdrf_file = Path(sdrf_file)
        self.num_files = num_files
        self.keep_files = keep_files
        self.file_types = file_types
        self.tolerance_unit = tolerance_unit
        self.skip_ptm = skip_ptm

        # Set up output paths
        stem = self.sdrf_file.stem.replace(".sdrf", "")
        parent = self.sdrf_file.parent

        self.output_sdrf = Path(output_sdrf) if output_sdrf else parent / f"{stem}_refined.sdrf.tsv"
        self.report_file = Path(report_file) if report_file else parent / f"{stem}_report.txt"

        # Set up working directory (will be cleaned up automatically)
        self.temp_manager = TempManager(
            work_dir=Path(work_dir) if work_dir else None,
            keep_files=keep_files,
        )

    def run(self) -> RefinementResult:
        """
        Execute the full refinement workflow.

        Returns:
            RefinementResult with outcome and output paths
        """
        try:
            with self.temp_manager:
                return self._run_workflow()
        except Exception as e:
            import traceback
            logger.error(f"Refinement failed: {e}")
            logger.error(traceback.format_exc())
            return RefinementResult(success=False, error_message=str(e))

    def _run_workflow(self) -> RefinementResult:
        """Execute workflow steps with incremental file processing.

        Files are downloaded, converted, and analyzed one at a time.
        By default, raw and mzML files are removed after each file is
        processed to minimise disk usage.  Use ``--keep-files`` to
        retain them (useful for resuming after pipeline changes).
        """
        # Step 1: Validate input SDRF
        logger.info("Step 1: Validating input SDRF...")
        sdrf_reader = SDRFReader(self.sdrf_file)

        validation_errors = sdrf_reader.validate()
        if validation_errors:
            logger.warning(f"SDRF validation found {len(validation_errors)} issues")
            for err in validation_errors[:5]:
                logger.warning(f"  - {err}")

        # Read SDRF and extract parameters
        sdrf_reader.read()
        sdrf_params = sdrf_reader.get_all_parameters()
        logger.info(f"Original SDRF parameters: {list(sdrf_params.keys())}")

        # Get data files from SDRF
        data_files = sdrf_reader.get_data_files()
        if not data_files:
            return RefinementResult(
                success=False,
                error_message="No data files found in SDRF comment[data file] column",
            )
        logger.info(f"Found {len(data_files)} data files in SDRF")

        # --- Resolve the list of filenames to process ---
        downloader = create_downloader(
            self.accession, self.temp_manager.download_dir
        )
        filenames_to_process = downloader.resolve_filenames(
            data_files, num_files=self.num_files, file_types=self.file_types
        )

        if not filenames_to_process:
            return RefinementResult(
                success=False,
                error_message="No raw files matched the requested criteria",
            )

        total = len(filenames_to_process)
        if self.num_files is None:
            logger.info(f"Step 2-4: Processing all {total} files incrementally...")
        else:
            logger.info(f"Step 2-4: Processing {total} files incrementally...")

        # --- Incremental per-file loop: download → convert → analyze → cleanup ---
        converter = RawFileConverter(self.temp_manager.mzml_dir)
        analyzer = MSAnalyzer(verbose=1)
        tolerance_estimator = ToleranceEstimator(use_gaussian=True, use_crux=True)
        ptm_detector = None if self.skip_ptm else PTMDetector()

        analysis_results: List[AnalysisResult] = []
        ptm_results: List[PTMDetectionResult] = []
        tolerance_result = None
        files_processed = 0

        for idx, filename in enumerate(filenames_to_process, 1):
            logger.info(f"[{idx}/{total}] Processing {filename}...")

            # -- Download --
            raw_file = downloader.download_file_by_name(filename)
            if not raw_file:
                logger.warning(f"[{idx}/{total}] Failed to download {filename}, skipping")
                continue

            # -- Convert --
            if raw_file.suffix.lower() == ".mzml":
                mzml_file = raw_file
                logger.info(f"[{idx}/{total}] Already mzML, skipping conversion")
            else:
                mzml_file = converter.convert(raw_file)
                if not mzml_file:
                    logger.warning(f"[{idx}/{total}] Conversion failed for {filename}, skipping")
                    self.temp_manager.remove_file(raw_file)
                    continue

            # -- Analyze --
            result = analyzer.analyze(mzml_file)
            analysis_results.append(result)
            files_processed += 1
            logger.info(
                f"[{idx}/{total}] Analyzed: instrument={result.instrument_model.get('name', 'unknown')}, "
                f"fragmentation={result.fragmentation_type}, spectra={result.n_spectra}"
            )

            # -- Tolerance estimation on first file only --
            if tolerance_result is None:
                logger.info(f"[{idx}/{total}] Estimating tolerances empirically (first file)...")
                tolerance_result = tolerance_estimator.estimate([mzml_file])
                logger.info(tolerance_result.format_summary(tolerance_result.reference_mz))

            # -- PTM detection --
            if ptm_detector is not None:
                logger.info(f"[{idx}/{total}] Detecting PTMs from spectra...")
                ptm_result = ptm_detector.detect(mzml_file)
                ptm_results.append(ptm_result)
                if ptm_result.hits:
                    hit_names = [h.name for h in ptm_result.hits]
                    logger.info(f"[{idx}/{total}] PTMs detected: {', '.join(hit_names)}")

            # -- Cleanup: remove raw + mzml to free disk space --
            if raw_file != mzml_file:
                self.temp_manager.remove_file(raw_file)
            self.temp_manager.remove_file(mzml_file)

        if not analysis_results:
            return RefinementResult(
                success=False,
                error_message="No files could be downloaded, converted, or analyzed",
            )

        logger.info(f"Processed {files_processed}/{total} files successfully")

        # Step 4: Aggregate analysis results
        logger.info("Step 4: Aggregating analysis results...")
        detected_params = analyzer.aggregate_results(analysis_results)
        logger.info(f"Detected parameters: {list(detected_params.keys())}")

        # Step 4b: Apply tolerance estimation to detected params
        reference_mz = tolerance_result.reference_mz if tolerance_result else 500.0

        # ── Precursor tolerance (always ppm, --tolerance-unit does NOT apply) ──
        if tolerance_result and tolerance_result.recommended_precursor_ppm is not None:
            detected_params["precursor_tolerance_ppm"] = tolerance_result.recommended_precursor_ppm
            detected_params["precursor_tolerance_source"] = "empirical"
            logger.info(f"Empirical precursor tolerance: {tolerance_result.recommended_precursor_ppm:.2f} ppm")
        else:
            detected_params["precursor_tolerance_source"] = "heuristic"
            logger.info("No empirical precursor tolerance available, using instrument heuristic")

        # ── Fragment tolerance (uses --tolerance-unit if set) ──
        if tolerance_result and tolerance_result.recommended_fragment_ppm is not None:
            detected_params["fragment_tolerance_ppm"] = tolerance_result.recommended_fragment_ppm
            detected_params["fragment_tolerance_da"] = tolerance_result.recommended_fragment_da
            detected_params["fragment_tolerance_source"] = "empirical"
            native = tolerance_result.fragment_native_unit or "ppm"
            detected_params["fragment_native_unit"] = native
            logger.info(
                f"Empirical fragment tolerance: "
                f"{tolerance_result.recommended_fragment_ppm:.2f} ppm "
                f"/ {tolerance_result.recommended_fragment_da:.4f} Da "
                f"(native unit from estimator: {native})"
            )
        else:
            detected_params["fragment_tolerance_source"] = "heuristic"
            logger.info("No empirical fragment tolerance available, using instrument heuristic")

        # Store reference m/z and tolerance estimation details
        detected_params["reference_mz"] = reference_mz
        detected_params["tolerance_estimation"] = tolerance_result

        # Step 4c: Aggregate PTM detection results
        if ptm_detector is not None and ptm_results:
            logger.info("Aggregating PTM detection results...")
            ptm_consensus = ptm_detector.aggregate_results(ptm_results)
            detected_params["ptm_detection"] = ptm_consensus
            n_ptms = len(ptm_consensus.get("hits", []))
            logger.info(f"PTM detection: {n_ptms} PTM types detected across {len(ptm_results)} files")

        # Step 5: Compare parameters
        logger.info("Step 5: Comparing detected vs SDRF parameters...")
        comparer = ParameterComparer(
            tolerance_unit=self.tolerance_unit,
            reference_mz=reference_mz,
        )
        comparisons = comparer.compare_all(sdrf_params, detected_params)
        refinements = comparer.get_refinements()

        logger.info(f"Found {len(refinements)} parameters that can be refined")

        # Step 6: Write refined SDRF (do this before report so we have the changes)
        logger.info("Step 6: Writing refined SDRF...")
        # Pass original columns to preserve duplicate column names
        original_columns = getattr(sdrf_reader, '_original_columns', list(sdrf_reader.df.columns))
        writer = SDRFWriter(sdrf_reader.df, original_columns=original_columns)
        num_applied = writer.apply_refinements(refinements)
        writer.write(self.output_sdrf)

        # Get the changes for the report
        changes = writer.get_changes()

        # Step 7: Generate report with changes
        logger.info("Step 7: Generating refinement report...")
        report_generator = ReportGenerator()
        report_generator.generate(
            pride_accession=self.pride_accession,
            sdrf_params=sdrf_params,
            detected_params=detected_params,
            comparisons=comparisons,
            analysis_results=analysis_results,
            output_path=self.report_file,
            changes=changes,
            ptm_results=ptm_results if ptm_results else None,
        )

        # Step 8: Validate refined SDRF
        logger.info("Step 8: Validating refined SDRF...")
        validation_errors = self._validate_refined_sdrf(
            self.output_sdrf, detected_params.get("acquisition_type", "DDA")
        )
        if validation_errors:
            logger.warning(f"Refined SDRF has {len(validation_errors)} validation issues:")
            for err in validation_errors[:10]:
                logger.warning(f"  - {err}")
            if len(validation_errors) > 10:
                logger.warning(f"  ... and {len(validation_errors) - 10} more")
        else:
            logger.info("Refined SDRF validation passed!")

        return RefinementResult(
            success=True,
            output_sdrf=self.output_sdrf,
            report_file=self.report_file,
            comparisons=comparisons,
            detected_params=detected_params,
            num_refinements=num_applied,
        )

    def _validate_refined_sdrf(
        self, sdrf_path: Path, acquisition_type: str
    ) -> List[str]:
        """
        Validate the refined SDRF using sdrf-pipelines.

        Uses the appropriate template based on acquisition type:
        - DDA: ms-proteomics + dda-acquisition
        - DIA: ms-proteomics + dia-acquisition

        Args:
            sdrf_path: Path to the refined SDRF file
            acquisition_type: Detected acquisition type (DDA/DIA)

        Returns:
            List of validation error messages
        """
        errors = []

        # Determine template based on acquisition type
        if acquisition_type == "DIA":
            template = "dia-acquisition"
        else:
            template = "dda-acquisition"

        logger.info(f"Validating with template: ms-proteomics + {template}")

        try:
            from sdrf_pipelines.sdrf.sdrf import read_sdrf

            sdrf_df = read_sdrf(str(sdrf_path))
            validation_errors = sdrf_df.validate_sdrf(template=template)

            for err in validation_errors:
                errors.append(str(err.message) if hasattr(err, "message") else str(err))

        except ImportError:
            logger.warning("sdrf-pipelines not available, skipping validation")
        except Exception as e:
            errors.append(f"Validation error: {str(e)}")

        return errors
