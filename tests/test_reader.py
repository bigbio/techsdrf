"""
Tests for SDRFReader - SDRF file reading and parsing.
"""

import pytest
from pathlib import Path

from sdrf_refiner.sdrf.reader import SDRFReader


FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestSDRFReader:
    """Tests for SDRFReader class."""

    def test_read_basic_sdrf(self):
        """Test reading a basic SDRF file."""
        reader = SDRFReader(FIXTURES_DIR / "sample_sdrf.tsv")
        df = reader.read()

        assert len(df) == 3
        assert "source name" in df.columns
        assert "comment[data file]" in df.columns

    def test_read_sdrf_with_metadata(self):
        """Test reading SDRF with metadata lines (#version, #project)."""
        reader = SDRFReader(FIXTURES_DIR / "sdrf_with_metadata.tsv")
        df = reader.read()

        assert len(df) == 1
        # Metadata should be parsed into _metadata dict
        assert reader._metadata.get("version") == "1.0"
        assert reader._metadata.get("project") == "PXD000001"

    def test_df_property_lazy_load(self):
        """Test that df property triggers read if not already loaded."""
        reader = SDRFReader(FIXTURES_DIR / "minimal_sdrf.tsv")
        assert reader._df is None

        # Accessing df should trigger read
        df = reader.df
        assert df is not None
        assert len(df) == 1

    def test_get_data_files(self):
        """Test extracting data file names from SDRF."""
        reader = SDRFReader(FIXTURES_DIR / "sample_sdrf.tsv")
        data_files = reader.get_data_files()

        assert len(data_files) == 3
        assert "sample1.raw" in data_files
        assert "sample2.raw" in data_files
        assert "sample3.raw" in data_files

    def test_get_parameter_for_file(self):
        """Test getting a specific parameter for a data file."""
        reader = SDRFReader(FIXTURES_DIR / "sample_sdrf.tsv")
        reader.read()

        instrument = reader.get_parameter_for_file(
            "sample1.raw", "comment[instrument]"
        )
        assert "Q Exactive HF" in instrument

    def test_get_parameter_for_file_missing_column(self):
        """Test getting parameter when column doesn't exist."""
        reader = SDRFReader(FIXTURES_DIR / "minimal_sdrf.tsv")
        reader.read()

        result = reader.get_parameter_for_file("sample1.raw", "comment[nonexistent]")
        assert result is None

    def test_get_all_parameters(self):
        """Test extracting all parameters from SDRF."""
        reader = SDRFReader(FIXTURES_DIR / "sample_sdrf.tsv")
        reader.read()
        params = reader.get_all_parameters()

        assert "instrument" in params
        assert params["instrument"]["name"] == "Q Exactive HF"
        assert params["instrument"]["accession"] == "MS:1002523"

        assert "dissociation" in params
        assert params["dissociation"]["name"] == "HCD"

        assert "precursor_tolerance" in params
        assert params["precursor_tolerance"]["value"] == 10.0
        assert params["precursor_tolerance"]["unit"] == "ppm"

        assert "fragment_tolerance" in params
        assert params["fragment_tolerance"]["value"] == 0.02
        assert params["fragment_tolerance"]["unit"] == "Da"

        assert "ms_min_charge" in params
        assert params["ms_min_charge"] == 2

        assert "ms_max_charge" in params
        assert params["ms_max_charge"] == 6

    def test_get_all_parameters_empty_sdrf(self):
        """Test getting parameters from minimal SDRF without parameter columns."""
        reader = SDRFReader(FIXTURES_DIR / "minimal_sdrf.tsv")
        reader.read()
        params = reader.get_all_parameters()

        # Should return empty dict or dict without those keys
        assert "instrument" not in params
        assert "dissociation" not in params


class TestOntologyValueParsing:
    """Tests for parsing ontology-formatted values."""

    def test_parse_full_ontology_value(self):
        """Test parsing NT=...;AC=... format."""
        result = SDRFReader.parse_ontology_value("NT=Q Exactive HF;AC=MS:1002523")

        assert result["name"] == "Q Exactive HF"
        assert result["accession"] == "MS:1002523"
        assert result["raw"] == "NT=Q Exactive HF;AC=MS:1002523"

    def test_parse_ontology_name_only(self):
        """Test parsing when only name is provided."""
        result = SDRFReader.parse_ontology_value("NT=Orbitrap Fusion")

        assert result["name"] == "Orbitrap Fusion"
        assert result["accession"] is None

    def test_parse_raw_value_no_format(self):
        """Test parsing plain value without NT/AC format."""
        result = SDRFReader.parse_ontology_value("Q Exactive")

        assert result["name"] == "Q Exactive"
        assert result["accession"] is None

    def test_parse_empty_value(self):
        """Test parsing empty or None value."""
        result = SDRFReader.parse_ontology_value("")
        assert result["name"] is None
        assert result["accession"] is None

        result = SDRFReader.parse_ontology_value(None)
        assert result["name"] is None
        assert result["accession"] is None


class TestToleranceValueParsing:
    """Tests for parsing tolerance values."""

    def test_parse_ppm_tolerance(self):
        """Test parsing ppm tolerance."""
        result = SDRFReader.parse_tolerance_value("10 ppm")

        assert result["value"] == 10.0
        assert result["unit"] == "ppm"

    def test_parse_da_tolerance(self):
        """Test parsing Da tolerance."""
        result = SDRFReader.parse_tolerance_value("0.02 Da")

        assert result["value"] == 0.02
        assert result["unit"] == "Da"

    def test_parse_tolerance_no_space(self):
        """Test parsing tolerance without space."""
        result = SDRFReader.parse_tolerance_value("20ppm")

        assert result["value"] == 20.0
        assert result["unit"] == "ppm"

    def test_parse_empty_tolerance(self):
        """Test parsing empty tolerance."""
        result = SDRFReader.parse_tolerance_value("")

        assert result["value"] is None
        assert result["unit"] is None

    def test_parse_tolerance_decimal(self):
        """Test parsing decimal tolerance values."""
        result = SDRFReader.parse_tolerance_value("0.5 Da")

        assert result["value"] == 0.5
        assert result["unit"] == "Da"
