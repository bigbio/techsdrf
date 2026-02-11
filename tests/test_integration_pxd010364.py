"""
Integration test for PXD010364 - Erythrocyte proteomics (TMT, Q Exactive HF).

Uses the SDRF from multiomics-configs. Marked as slow because it downloads
real data from PRIDE, converts, and analyzes.
"""

import logging
import shutil
import tempfile
from pathlib import Path

import pytest

# Path to the PXD010364 SDRF in multiomics-configs
MULTIOMICS_SDRF = (
    Path(__file__).parent.parent.parent
    / "multiomics-configs"
    / "projects"
    / "blood-proteomes"
    / "Erythrocytes"
    / "PXD010364"
    / "PXD010364.sdrf.tsv"
)


class TestPXD010364Integration:
    """End-to-end integration tests for PXD010364."""

    @pytest.fixture(autouse=True)
    def setup_work_dir(self, tmp_path):
        """Create a temporary work directory for each test."""
        self.work_dir = tmp_path / "pxd010364_test"
        self.work_dir.mkdir()

    def test_sdrf_file_exists(self):
        """Verify the PXD010364 SDRF file is available."""
        assert MULTIOMICS_SDRF.exists(), (
            f"PXD010364 SDRF not found at {MULTIOMICS_SDRF}. "
            "Make sure multiomics-configs is checked out."
        )

    def test_sdrf_reader_parses_pxd010364(self):
        """Test that SDRFReader can parse the PXD010364 SDRF."""
        if not MULTIOMICS_SDRF.exists():
            pytest.skip("multiomics-configs not available")

        from sdrf_refiner.sdrf.reader import SDRFReader

        reader = SDRFReader(MULTIOMICS_SDRF)
        df = reader.read()

        # PXD010364 has 20 rows (10 per raw file, 2 raw files in default SDRF)
        assert len(df) >= 20
        assert "comment[data file]" in df.columns
        assert "comment[instrument]" in df.columns

        # Check data files are .raw
        data_files = reader.get_data_files()
        assert len(data_files) > 0
        assert all(f.endswith(".raw") for f in data_files)

        # Check instrument is Q Exactive HF
        params = reader.get_all_parameters()
        instrument = params.get("instrument", {})
        assert "Q Exactive HF" in str(instrument) or "MS:1002523" in str(instrument)

    def test_sdrf_parameters_extraction(self):
        """Verify key parameters are correctly extracted from PXD010364 SDRF."""
        if not MULTIOMICS_SDRF.exists():
            pytest.skip("multiomics-configs not available")

        from sdrf_refiner.sdrf.reader import SDRFReader

        reader = SDRFReader(MULTIOMICS_SDRF)
        reader.read()
        params = reader.get_all_parameters()

        # Instrument: Q Exactive HF (MS:1002523)
        assert "instrument" in params
        assert "MS:1002523" in str(params["instrument"])

        # Dissociation: HCD (MS:1000422)
        assert "dissociation" in params
        assert "MS:1000422" in str(params["dissociation"])

        # Precursor tolerance: 5 ppm
        assert "precursor_tolerance" in params

        # Fragment tolerance: 0.2 Da
        assert "fragment_tolerance" in params

    def test_downloader_resolves_filenames(self):
        """Test that PrideDownloader resolves filenames from SDRF."""
        if not MULTIOMICS_SDRF.exists():
            pytest.skip("multiomics-configs not available")

        from sdrf_refiner.pride.downloader import PrideDownloader
        from sdrf_refiner.sdrf.reader import SDRFReader

        reader = SDRFReader(MULTIOMICS_SDRF)
        reader.read()
        data_files = reader.get_data_files()

        downloader = PrideDownloader("PXD010364", self.work_dir / "downloads")

        # Test resolve_filenames with limit
        filenames = downloader.resolve_filenames(data_files, num_files=2)
        assert len(filenames) == 2
        assert all(f.endswith(".raw") for f in filenames)

        # Test resolve_filenames without limit (should get all unique)
        all_filenames = downloader.resolve_filenames(data_files)
        # PXD010364 SDRF references multiple .raw files
        assert len(all_filenames) >= 2
        # No duplicates
        assert len(all_filenames) == len(set(all_filenames))

    @pytest.mark.slow
    def test_refine_pxd010364_single_file(self):
        """
        Full end-to-end refine for PXD010364 with 1 file.

        Downloads, converts, analyzes, and produces refined SDRF.
        Verifies incremental processing (files cleaned up after each).
        """
        if not MULTIOMICS_SDRF.exists():
            pytest.skip("multiomics-configs not available")

        # Check ThermoRawFileParser is available
        if not shutil.which("ThermoRawFileParser") and not shutil.which("ThermoRawFileParser.sh"):
            pytest.skip("ThermoRawFileParser not available")

        from sdrf_refiner.core.refiner import SDRFRefiner

        output_sdrf = self.work_dir / "PXD010364_refined.sdrf.tsv"
        report_file = self.work_dir / "PXD010364_report.txt"

        refiner = SDRFRefiner(
            pride_accession="PXD010364",
            sdrf_file=MULTIOMICS_SDRF,
            num_files=1,
            output_sdrf=output_sdrf,
            report_file=report_file,
            work_dir=self.work_dir,
            keep_files=False,  # Default: clean up after each file
        )

        result = refiner.run()

        assert result.success, f"Refinement failed: {result.error_message}"
        assert output_sdrf.exists(), "Refined SDRF not created"
        assert report_file.exists(), "Report not created"

        # Verify files were cleaned up (incremental processing)
        downloads_dir = self.work_dir / "downloads"
        mzml_dir = self.work_dir / "mzml"
        raw_files = list(downloads_dir.glob("*.raw"))
        mzml_files = list(mzml_dir.glob("*.mzML"))
        assert len(raw_files) == 0, f"Raw files not cleaned up: {raw_files}"
        assert len(mzml_files) == 0, f"mzML files not cleaned up: {mzml_files}"

    @pytest.mark.slow
    def test_refine_pxd010364_keep_files(self):
        """
        Full end-to-end refine with --keep-files.

        Verifies files are retained after processing.
        """
        if not MULTIOMICS_SDRF.exists():
            pytest.skip("multiomics-configs not available")

        if not shutil.which("ThermoRawFileParser") and not shutil.which("ThermoRawFileParser.sh"):
            pytest.skip("ThermoRawFileParser not available")

        from sdrf_refiner.core.refiner import SDRFRefiner

        output_sdrf = self.work_dir / "PXD010364_refined.sdrf.tsv"
        report_file = self.work_dir / "PXD010364_report.txt"

        refiner = SDRFRefiner(
            pride_accession="PXD010364",
            sdrf_file=MULTIOMICS_SDRF,
            num_files=1,
            output_sdrf=output_sdrf,
            report_file=report_file,
            work_dir=self.work_dir,
            keep_files=True,  # Keep files for resume
        )

        result = refiner.run()

        assert result.success, f"Refinement failed: {result.error_message}"

        # Verify files were kept
        downloads_dir = self.work_dir / "downloads"
        mzml_dir = self.work_dir / "mzml"
        raw_files = list(downloads_dir.glob("*.raw"))
        mzml_files = list(mzml_dir.glob("*.mzML"))
        assert len(raw_files) >= 1, "Raw files should be kept with --keep-files"
        assert len(mzml_files) >= 1, "mzML files should be kept with --keep-files"
