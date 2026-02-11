"""
Tests for SDRFWriter - SDRF file writing and refinement application.
"""

import pytest
import pandas as pd
from io import StringIO
from pathlib import Path
import tempfile

from sdrf_refiner.sdrf.writer import SDRFWriter
from sdrf_refiner.sdrf.reader import SDRFReader
from sdrf_refiner.sdrf.comparer import ParameterComparison, ComparisonStatus


FIXTURES_DIR = Path(__file__).parent / "fixtures"


class TestSDRFWriterInit:
    """Tests for SDRFWriter initialization."""

    def test_init_copies_dataframe(self):
        """Test that initialization makes a copy of the DataFrame."""
        df = pd.DataFrame({"col1": ["a", "b"], "col2": [1, 2]})
        writer = SDRFWriter(df)

        # Modify original
        df["col1"] = ["x", "y"]

        # Writer's copy should be unchanged
        assert writer.df["col1"].tolist() == ["a", "b"]

    def test_init_with_custom_columns(self):
        """Test initialization with custom original column names."""
        df = pd.DataFrame({"col1": ["a"], "col2": ["b"]})
        original_cols = ["col1", "col2", "col2"]  # Duplicate column name

        writer = SDRFWriter(df, original_columns=original_cols)

        assert writer._original_columns == original_cols


class TestApplyRefinements:
    """Tests for applying refinements to SDRF."""

    def test_apply_mismatch_refinement(self):
        """Test applying a mismatch refinement."""
        df = pd.DataFrame({
            "source name": ["sample1"],
            "comment[instrument]": ["NT=Old Instrument;AC=MS:0000001"],
        })

        writer = SDRFWriter(df)
        refinements = [
            ParameterComparison(
                parameter_name="instrument",
                sdrf_value={"name": "Old Instrument"},
                detected_value={"name": "Q Exactive HF", "accession": "MS:1002523"},
                status=ComparisonStatus.MISMATCH,
                recommendation="Update instrument",
            ),
        ]

        applied = writer.apply_refinements(refinements)

        assert applied == 1
        assert "Q Exactive HF" in writer.df["comment[instrument]"].iloc[0]

    def test_apply_missing_sdrf_refinement(self):
        """Test applying refinement when SDRF column is missing."""
        df = pd.DataFrame({
            "source name": ["sample1"],
        })

        writer = SDRFWriter(df)
        refinements = [
            ParameterComparison(
                parameter_name="instrument",
                sdrf_value=None,
                detected_value={"name": "Q Exactive HF", "accession": "MS:1002523"},
                status=ComparisonStatus.MISSING_SDRF,
                recommendation="Add instrument",
            ),
        ]

        applied = writer.apply_refinements(refinements)

        assert applied == 1
        assert "comment[instrument]" in writer.df.columns
        assert "Q Exactive HF" in writer.df["comment[instrument]"].iloc[0]

    def test_apply_tolerance_refinement(self):
        """Test applying tolerance refinement."""
        df = pd.DataFrame({
            "source name": ["sample1"],
            "comment[precursor mass tolerance]": ["20 ppm"],
        })

        writer = SDRFWriter(df)
        refinements = [
            ParameterComparison(
                parameter_name="precursor_mass_tolerance",
                sdrf_value={"value": 20, "unit": "ppm"},
                detected_value={"value": 10.7, "unit": "ppm"},
                status=ComparisonStatus.IMPROVED,
                recommendation="Update tolerance",
            ),
        ]

        applied = writer.apply_refinements(refinements)

        assert applied == 1
        # Should round to integer for ppm (10.7 rounds to 11)
        assert "11 ppm" in writer.df["comment[precursor mass tolerance]"].iloc[0]

    def test_apply_charge_refinement(self):
        """Test applying charge state refinement."""
        df = pd.DataFrame({
            "source name": ["sample1"],
            "comment[MS min charge]": ["1"],
        })

        writer = SDRFWriter(df)
        refinements = [
            ParameterComparison(
                parameter_name="ms_min_charge",
                sdrf_value=1,
                detected_value=2,
                status=ComparisonStatus.MISMATCH,
                recommendation="Update charge",
            ),
        ]

        applied = writer.apply_refinements(refinements)

        assert applied == 1
        assert writer.df["comment[MS min charge]"].iloc[0] == "2"

    def test_skip_match_refinement(self):
        """Test that MATCH status is skipped."""
        df = pd.DataFrame({
            "source name": ["sample1"],
            "comment[instrument]": ["NT=Q Exactive HF;AC=MS:1002523"],
        })

        writer = SDRFWriter(df)
        refinements = [
            ParameterComparison(
                parameter_name="instrument",
                sdrf_value={"name": "Q Exactive HF"},
                detected_value={"name": "Q Exactive HF"},
                status=ComparisonStatus.MATCH,
            ),
        ]

        applied = writer.apply_refinements(refinements)

        assert applied == 0

    def test_track_changes(self):
        """Test that changes are tracked for reporting."""
        df = pd.DataFrame({
            "source name": ["sample1"],
            "comment[instrument]": ["NT=Old;AC=MS:0000001"],
        })

        writer = SDRFWriter(df)
        refinements = [
            ParameterComparison(
                parameter_name="instrument",
                sdrf_value={"name": "Old"},
                detected_value={"name": "New", "accession": "MS:1002523"},
                status=ComparisonStatus.MISMATCH,
                recommendation="Update instrument",
            ),
        ]

        writer.apply_refinements(refinements)
        changes = writer.get_changes()

        assert len(changes) == 1
        assert changes[0]["column"] == "comment[instrument]"
        assert "Old" in changes[0]["old_value"]
        assert "New" in changes[0]["new_value"]


class TestWrite:
    """Tests for writing SDRF files."""

    def test_write_basic_sdrf(self):
        """Test writing a basic SDRF file."""
        df = pd.DataFrame({
            "source name": ["sample1", "sample2"],
            "comment[data file]": ["file1.raw", "file2.raw"],
        })

        writer = SDRFWriter(df)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.sdrf.tsv"
            writer.write(output_path)

            assert output_path.exists()

            # Read back and verify
            with open(output_path, "r") as f:
                lines = f.readlines()

            assert len(lines) == 3  # header + 2 data rows
            assert "source name" in lines[0]
            assert "sample1" in lines[1]
            assert "sample2" in lines[2]

    def test_write_creates_parent_dirs(self):
        """Test that write creates parent directories if needed."""
        df = pd.DataFrame({"col1": ["a"]})
        writer = SDRFWriter(df)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "subdir" / "nested" / "output.sdrf.tsv"
            writer.write(output_path)

            assert output_path.exists()


class TestFormatValue:
    """Tests for _format_value static method."""

    def test_format_instrument_value(self):
        """Test formatting instrument value."""
        comparison = ParameterComparison(
            parameter_name="instrument",
            sdrf_value=None,
            detected_value={"name": "Q Exactive HF", "accession": "MS:1002523"},
            status=ComparisonStatus.MISSING_SDRF,
        )

        result = SDRFWriter._format_value(comparison)

        assert result == "AC=MS:1002523;NT=Q Exactive HF"

    def test_format_instrument_name_only(self):
        """Test formatting instrument with name only."""
        comparison = ParameterComparison(
            parameter_name="instrument",
            sdrf_value=None,
            detected_value={"name": "Q Exactive HF"},
            status=ComparisonStatus.MISSING_SDRF,
        )

        result = SDRFWriter._format_value(comparison)

        assert result == "Q Exactive HF"

    def test_format_tolerance_ppm(self):
        """Test formatting ppm tolerance."""
        comparison = ParameterComparison(
            parameter_name="precursor_mass_tolerance",
            sdrf_value=None,
            detected_value={"value": 10.7, "unit": "ppm"},
            status=ComparisonStatus.MISSING_SDRF,
        )

        result = SDRFWriter._format_value(comparison)

        assert result == "11 ppm"  # Rounded to integer

    def test_format_tolerance_da(self):
        """Test formatting Da tolerance."""
        comparison = ParameterComparison(
            parameter_name="fragment_mass_tolerance",
            sdrf_value=None,
            detected_value={"value": 0.02, "unit": "Da"},
            status=ComparisonStatus.MISSING_SDRF,
        )

        result = SDRFWriter._format_value(comparison)

        assert "0.02" in result
        assert "Da" in result

    def test_format_charge_value(self):
        """Test formatting charge value."""
        comparison = ParameterComparison(
            parameter_name="ms_min_charge",
            sdrf_value=None,
            detected_value=2,
            status=ComparisonStatus.MISSING_SDRF,
        )

        result = SDRFWriter._format_value(comparison)

        assert result == "2"

    def test_format_dissociation_method(self):
        """Test formatting dissociation method."""
        comparison = ParameterComparison(
            parameter_name="dissociation_method",
            sdrf_value=None,
            detected_value={"name": "HCD", "accession": "MS:1000422"},
            status=ComparisonStatus.MISSING_SDRF,
        )

        result = SDRFWriter._format_value(comparison)

        assert result == "AC=MS:1000422;NT=HCD"

    def test_format_empty_detected_value(self):
        """Test formatting with empty detected value."""
        comparison = ParameterComparison(
            parameter_name="instrument",
            sdrf_value=None,
            detected_value=None,
            status=ComparisonStatus.MISSING_DETECTED,
        )

        result = SDRFWriter._format_value(comparison)

        assert result is None


class TestParamToColumn:
    """Tests for _param_to_column static method."""

    def test_instrument_mapping(self):
        """Test instrument parameter mapping."""
        assert SDRFWriter._param_to_column("instrument") == "comment[instrument]"

    def test_tolerance_mapping(self):
        """Test tolerance parameter mappings."""
        assert SDRFWriter._param_to_column("precursor_mass_tolerance") == "comment[precursor mass tolerance]"
        assert SDRFWriter._param_to_column("fragment_mass_tolerance") == "comment[fragment mass tolerance]"

    def test_charge_mapping(self):
        """Test charge parameter mappings."""
        assert SDRFWriter._param_to_column("ms_min_charge") == "comment[MS min charge]"
        assert SDRFWriter._param_to_column("ms_max_charge") == "comment[MS max charge]"

    def test_unknown_param(self):
        """Test unknown parameter returns None."""
        assert SDRFWriter._param_to_column("unknown_param") is None


class TestScanWindowParams:
    """Tests for scan window and isolation width writing."""

    def test_param_to_column_scan_windows(self):
        assert SDRFWriter._param_to_column("scan_window_lower") == "comment[scan window lower limit]"
        assert SDRFWriter._param_to_column("scan_window_upper") == "comment[scan window upper limit]"
        assert SDRFWriter._param_to_column("isolation_window_width") == "comment[isolation window width]"

    def test_format_scan_window_integer(self):
        comp = ParameterComparison(
            parameter_name="scan_window_lower",
            sdrf_value=None,
            detected_value=350.0,
            status=ComparisonStatus.MISSING_SDRF,
        )
        assert SDRFWriter._format_value(comp) == "350"

    def test_format_scan_window_decimal(self):
        comp = ParameterComparison(
            parameter_name="scan_window_upper",
            sdrf_value=None,
            detected_value=1600.5,
            status=ComparisonStatus.MISSING_SDRF,
        )
        assert SDRFWriter._format_value(comp) == "1600.5"

    def test_format_isolation_window(self):
        comp = ParameterComparison(
            parameter_name="isolation_window_width",
            sdrf_value=None,
            detected_value=2.0,
            status=ComparisonStatus.MISSING_SDRF,
        )
        assert SDRFWriter._format_value(comp) == "2"

    def test_apply_scan_window_refinements(self):
        """Test that scan window columns are actually written to output."""
        df = pd.DataFrame({
            "source name": ["sample1"],
            "comment[data file]": ["sample1.raw"],
        })
        writer = SDRFWriter(df)

        refinements = [
            ParameterComparison(
                parameter_name="scan_window_lower",
                sdrf_value=None,
                detected_value=350.0,
                status=ComparisonStatus.MISSING_SDRF,
            ),
            ParameterComparison(
                parameter_name="scan_window_upper",
                sdrf_value=None,
                detected_value=1600.0,
                status=ComparisonStatus.MISSING_SDRF,
            ),
            ParameterComparison(
                parameter_name="isolation_window_width",
                sdrf_value=None,
                detected_value=25.0,
                status=ComparisonStatus.MISSING_SDRF,
            ),
        ]
        applied = writer.apply_refinements(refinements)
        assert applied == 3

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.sdrf.tsv"
            writer.write(output_path)

            with open(output_path, "r") as f:
                lines = f.readlines()

            header = lines[0].strip().split("\t")
            row1 = lines[1].strip().split("\t")

            assert "comment[scan window lower limit]" in header
            assert "comment[scan window upper limit]" in header
            assert "comment[isolation window width]" in header

            sw_lower_pos = header.index("comment[scan window lower limit]")
            sw_upper_pos = header.index("comment[scan window upper limit]")
            iw_pos = header.index("comment[isolation window width]")

            assert row1[sw_lower_pos] == "350"
            assert row1[sw_upper_pos] == "1600"
            assert row1[iw_pos] == "25"


class TestDuplicateColumns:
    """Tests for proper handling of duplicate column names (e.g. modification parameters)."""

    def test_write_preserves_duplicate_column_values(self):
        """Test that duplicate columns each keep their own distinct values."""
        # Simulate what pandas does when reading duplicate column names:
        # the first stays as-is, the second gets a .1 suffix
        df = pd.DataFrame({
            "source name": ["sample1", "sample2"],
            "comment[modification parameters]": [
                "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed",
                "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed",
            ],
            "comment[modification parameters].1": [
                "NT=Oxidation;MT=Variable;TA=M;AC=UNIMOD:35",
                "NT=Oxidation;MT=Variable;TA=M;AC=UNIMOD:35",
            ],
        })

        # Original columns as they appear in the real TSV (no .1 suffix)
        original_columns = [
            "source name",
            "comment[modification parameters]",
            "comment[modification parameters]",
        ]

        writer = SDRFWriter(df, original_columns=original_columns)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.sdrf.tsv"
            writer.write(output_path)

            with open(output_path, "r") as f:
                lines = f.readlines()

            # Header should have the duplicate name twice
            header_cols = lines[0].strip().split("\t")
            assert header_cols.count("comment[modification parameters]") == 2

            # Data: first column should be Carbamidomethyl, second should be Oxidation
            row1_cols = lines[1].strip().split("\t")
            assert "Carbamidomethyl" in row1_cols[1]
            assert "Oxidation" in row1_cols[2]

    def test_write_preserves_duplicates_after_refinement(self):
        """Test that applying refinements to other columns does not corrupt duplicate columns."""
        df = pd.DataFrame({
            "source name": ["sample1"],
            "comment[instrument]": ["NT=Q Exactive;AC=MS:1001911"],
            "comment[modification parameters]": [
                "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed",
            ],
            "comment[modification parameters].1": [
                "NT=Oxidation;MT=Variable;TA=M;AC=UNIMOD:35",
            ],
            "comment[precursor mass tolerance]": ["20 ppm"],
        })

        original_columns = [
            "source name",
            "comment[instrument]",
            "comment[modification parameters]",
            "comment[modification parameters]",
            "comment[precursor mass tolerance]",
        ]

        writer = SDRFWriter(df, original_columns=original_columns)

        # Apply tolerance refinement (should not touch modification parameters)
        refinements = [
            ParameterComparison(
                parameter_name="precursor_mass_tolerance",
                sdrf_value={"value": 20, "unit": "ppm"},
                detected_value={"value": 10.0, "unit": "ppm"},
                status=ComparisonStatus.IMPROVED,
                recommendation="Update tolerance",
            ),
        ]
        writer.apply_refinements(refinements)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.sdrf.tsv"
            writer.write(output_path)

            with open(output_path, "r") as f:
                lines = f.readlines()

            row1_cols = lines[1].strip().split("\t")
            # Modification columns must be preserved independently
            assert "Carbamidomethyl" in row1_cols[2]
            assert "Oxidation" in row1_cols[3]
            # Tolerance should be updated
            assert "10 ppm" in row1_cols[4]

    def test_roundtrip_reader_writer_preserves_duplicates(self):
        """End-to-end: read SDRF with duplicate columns, apply refinement, write back."""
        fixture_path = FIXTURES_DIR / "sdrf_with_duplicate_columns.tsv"

        reader = SDRFReader(fixture_path)
        reader.read()

        original_columns = reader._original_columns
        writer = SDRFWriter(reader.df, original_columns=original_columns)

        # Apply a tolerance refinement
        refinements = [
            ParameterComparison(
                parameter_name="precursor_mass_tolerance",
                sdrf_value={"value": 20, "unit": "ppm"},
                detected_value={"value": 7.0, "unit": "ppm"},
                status=ComparisonStatus.IMPROVED,
                recommendation="Update tolerance",
            ),
        ]
        writer.apply_refinements(refinements)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.sdrf.tsv"
            writer.write(output_path)

            # Read back the output and verify
            with open(output_path, "r") as f:
                lines = f.readlines()

            header = lines[0].strip().split("\t")
            row1 = lines[1].strip().split("\t")
            row2 = lines[2].strip().split("\t")

            # Header must have duplicate column names
            mod_param_positions = [
                i for i, h in enumerate(header) if h == "comment[modification parameters]"
            ]
            assert len(mod_param_positions) == 2

            # First mod param column: Carbamidomethyl
            assert "Carbamidomethyl" in row1[mod_param_positions[0]]
            assert "Carbamidomethyl" in row2[mod_param_positions[0]]

            # Second mod param column: Oxidation
            assert "Oxidation" in row1[mod_param_positions[1]]
            assert "Oxidation" in row2[mod_param_positions[1]]

            # Tolerance should be updated
            tol_pos = header.index("comment[precursor mass tolerance]")
            assert "7 ppm" in row1[tol_pos]

    def test_adding_new_column_with_existing_duplicates(self):
        """Test adding a new column (e.g. dissociation method) when duplicates exist."""
        df = pd.DataFrame({
            "source name": ["sample1"],
            "comment[modification parameters]": [
                "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed",
            ],
            "comment[modification parameters].1": [
                "NT=Oxidation;MT=Variable;TA=M;AC=UNIMOD:35",
            ],
        })

        original_columns = [
            "source name",
            "comment[modification parameters]",
            "comment[modification parameters]",
        ]

        writer = SDRFWriter(df, original_columns=original_columns)

        # Add a new column via refinement
        refinements = [
            ParameterComparison(
                parameter_name="dissociation_method",
                sdrf_value=None,
                detected_value={"name": "HCD", "accession": "MS:1000422"},
                status=ComparisonStatus.MISSING_SDRF,
                recommendation="Add dissociation method",
            ),
        ]
        writer.apply_refinements(refinements)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.sdrf.tsv"
            writer.write(output_path)

            with open(output_path, "r") as f:
                lines = f.readlines()

            header = lines[0].strip().split("\t")
            row1 = lines[1].strip().split("\t")

            # Both modification columns should be intact
            mod_positions = [
                i for i, h in enumerate(header) if h == "comment[modification parameters]"
            ]
            assert len(mod_positions) == 2
            assert "Carbamidomethyl" in row1[mod_positions[0]]
            assert "Oxidation" in row1[mod_positions[1]]

            # New column should be present
            assert "comment[dissociation method]" in header
            diss_pos = header.index("comment[dissociation method]")
            assert "HCD" in row1[diss_pos]

    def test_no_original_columns_strips_pandas_suffixes(self):
        """Test that when original_columns is not provided, pandas .1/.2 suffixes are stripped."""
        df = pd.DataFrame({
            "source name": ["sample1"],
            "comment[modification parameters]": [
                "NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=Fixed",
            ],
            "comment[modification parameters].1": [
                "NT=Oxidation;MT=Variable;TA=M;AC=UNIMOD:35",
            ],
            "comment[precursor mass tolerance]": ["20 ppm"],
        })

        # Do NOT pass original_columns — writer must derive clean names
        writer = SDRFWriter(df)

        with tempfile.TemporaryDirectory() as tmpdir:
            output_path = Path(tmpdir) / "output.sdrf.tsv"
            writer.write(output_path)

            with open(output_path, "r") as f:
                lines = f.readlines()

            header = lines[0].strip().split("\t")
            row1 = lines[1].strip().split("\t")

            # Header must NOT contain .1 suffix
            assert "comment[modification parameters].1" not in header
            # Should have the duplicate name twice
            assert header.count("comment[modification parameters]") == 2

            # Values must be distinct
            mod_positions = [
                i for i, h in enumerate(header) if h == "comment[modification parameters]"
            ]
            assert "Carbamidomethyl" in row1[mod_positions[0]]
            assert "Oxidation" in row1[mod_positions[1]]

    def test_strip_pandas_suffixes_only_for_real_duplicates(self):
        """Test that .1 suffix is only stripped when the base name exists earlier."""
        # "score.1" is NOT a pandas duplicate — it's a legitimate column name
        df = pd.DataFrame({
            "source name": ["sample1"],
            "score.1": ["0.95"],
        })

        writer = SDRFWriter(df)

        # "score.1" should be kept as-is because "score" never appeared
        assert "score.1" in writer._original_columns

    def test_strip_pandas_suffixes_static(self):
        """Test _strip_pandas_suffixes directly."""
        # Standard duplicate case
        assert SDRFWriter._strip_pandas_suffixes(
            ["A", "B", "B.1", "C"]
        ) == ["A", "B", "B", "C"]

        # Multiple duplicates
        assert SDRFWriter._strip_pandas_suffixes(
            ["A", "A.1", "A.2"]
        ) == ["A", "A", "A"]

        # No duplicates — nothing stripped
        assert SDRFWriter._strip_pandas_suffixes(
            ["A", "B", "C"]
        ) == ["A", "B", "C"]

        # Legitimate .1 name (no base column)
        assert SDRFWriter._strip_pandas_suffixes(
            ["score.1", "name"]
        ) == ["score.1", "name"]

        # Mixed: one real duplicate, one legitimate .1
        assert SDRFWriter._strip_pandas_suffixes(
            ["mod", "mod.1", "score.1"]
        ) == ["mod", "mod", "score.1"]
