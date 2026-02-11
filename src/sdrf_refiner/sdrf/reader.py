"""
SDRFReader - Read and parse SDRF files using sdrf-pipelines.
"""

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from sdrf_refiner.config import SDRF_COLUMNS

logger = logging.getLogger(__name__)


class SDRFReader:
    """
    Reads and parses SDRF files, extracting parameter values
    for comparison with empirically detected values.
    """

    def __init__(self, sdrf_path: Path):
        """
        Initialize SDRF reader.

        Args:
            sdrf_path: Path to the SDRF file
        """
        self.sdrf_path = Path(sdrf_path)
        self._df: Optional[pd.DataFrame] = None
        self._metadata: Dict[str, str] = {}

    def read(self) -> pd.DataFrame:
        """
        Read the SDRF file and return as DataFrame.

        Returns:
            pandas DataFrame with SDRF contents
        """
        logger.info(f"Reading SDRF file: {self.sdrf_path}")

        # Read metadata lines (starting with #) and data
        metadata_lines = []
        data_lines = []

        with open(self.sdrf_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    metadata_lines.append(line.strip())
                else:
                    data_lines.append(line)

        # Parse metadata
        for meta_line in metadata_lines:
            if "=" in meta_line:
                key, value = meta_line.lstrip("#").split("=", 1)
                self._metadata[key.strip()] = value.strip()

        # Read data as DataFrame
        # IMPORTANT: SDRF allows duplicate column names (e.g., multiple comment[modification parameters])
        # We need to preserve these without pandas adding .1, .2 suffixes
        from io import StringIO

        # First, read header to get original column names
        header_line = data_lines[0] if data_lines else ""
        self._original_columns = header_line.strip().split("\t") if header_line else []

        # Read with pandas but preserve duplicate columns
        self._df = pd.read_csv(
            StringIO("".join(data_lines)), sep="\t", dtype=str, header=0
        ).fillna("")

        # Store the original column names for writing back
        # pandas may have renamed duplicates to .1, .2, etc.

        logger.info(f"Loaded SDRF with {len(self._df)} rows and {len(self._df.columns)} columns")
        return self._df

    @property
    def df(self) -> pd.DataFrame:
        """Get the DataFrame, reading the file if necessary."""
        if self._df is None:
            self.read()
        return self._df

    def validate(self) -> List[str]:
        """
        Validate the SDRF file using sdrf-pipelines.

        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []
        try:
            from sdrf_pipelines.sdrf.sdrf import read_sdrf

            sdrf_df = read_sdrf(str(self.sdrf_path))
            validation_errors = sdrf_df.validate_sdrf(skip_ontology=True)

            for err in validation_errors:
                errors.append(str(err.message) if hasattr(err, "message") else str(err))

        except ImportError:
            logger.warning("sdrf-pipelines not available, skipping validation")
        except Exception as e:
            errors.append(f"Validation error: {str(e)}")

        return errors

    def get_data_files(self) -> List[str]:
        """Get list of data files referenced in the SDRF."""
        data_file_col = SDRF_COLUMNS["data_file"]
        if data_file_col in self.df.columns:
            return self.df[data_file_col].unique().tolist()
        return []

    def get_parameter_for_file(self, data_file: str, column: str) -> Optional[str]:
        """
        Get a specific parameter value for a given data file.

        Args:
            data_file: The data file name
            column: The column name to retrieve

        Returns:
            The value or None if not found
        """
        if column not in self.df.columns:
            return None

        data_file_col = SDRF_COLUMNS["data_file"]
        if data_file_col not in self.df.columns:
            return None

        matching_rows = self.df[self.df[data_file_col] == data_file]
        if not matching_rows.empty:
            return matching_rows[column].iloc[0]
        return None

    def get_all_parameters(self) -> Dict[str, Any]:
        """
        Extract all relevant parameters from the SDRF.

        Returns:
            Dictionary with parsed parameter values
        """
        params = {}

        # Get parameters from first row (assuming consistent across file)
        if len(self.df) == 0:
            return params

        row = self.df.iloc[0]

        # Instrument
        inst_col = SDRF_COLUMNS["instrument"]
        if inst_col in self.df.columns and row[inst_col]:
            params["instrument"] = self.parse_ontology_value(row[inst_col])

        # Dissociation method (generic)
        diss_col = SDRF_COLUMNS["dissociation"]
        if diss_col in self.df.columns and row[diss_col]:
            params["dissociation"] = self.parse_ontology_value(row[diss_col])

        # MS2 dissociation method (specific)
        ms2_diss_col = SDRF_COLUMNS["ms2_dissociation"]
        if ms2_diss_col in self.df.columns and row[ms2_diss_col]:
            params["ms2_dissociation"] = self.parse_ontology_value(row[ms2_diss_col])

        # MS3 dissociation method (for SPS-MS3 experiments)
        ms3_diss_col = SDRF_COLUMNS["ms3_dissociation"]
        if ms3_diss_col in self.df.columns and row[ms3_diss_col]:
            params["ms3_dissociation"] = self.parse_ontology_value(row[ms3_diss_col])

        # Precursor tolerance
        prec_col = SDRF_COLUMNS["precursor_tolerance"]
        if prec_col in self.df.columns and row[prec_col]:
            params["precursor_tolerance"] = self.parse_tolerance_value(row[prec_col])

        # Fragment tolerance
        frag_col = SDRF_COLUMNS["fragment_tolerance"]
        if frag_col in self.df.columns and row[frag_col]:
            params["fragment_tolerance"] = self.parse_tolerance_value(row[frag_col])

        # Acquisition method
        acq_col = SDRF_COLUMNS["acquisition_method"]
        if acq_col in self.df.columns and row[acq_col]:
            params["acquisition_method"] = self.parse_ontology_value(row[acq_col])

        # MS min charge
        min_charge_col = SDRF_COLUMNS["ms_min_charge"]
        if min_charge_col in self.df.columns and row[min_charge_col]:
            try:
                params["ms_min_charge"] = int(row[min_charge_col])
            except (ValueError, TypeError):
                pass

        # MS max charge
        max_charge_col = SDRF_COLUMNS["ms_max_charge"]
        if max_charge_col in self.df.columns and row[max_charge_col]:
            try:
                params["ms_max_charge"] = int(row[max_charge_col])
            except (ValueError, TypeError):
                pass

        # Collision energy (format: "30 NCE", "27 eV", "25 NCE;27 NCE;30 NCE")
        ce_col = SDRF_COLUMNS["collision_energy"]
        if ce_col in self.df.columns and row[ce_col]:
            params["collision_energy"] = str(row[ce_col]).strip()

        # MS2 mass analyzer
        ma_col = SDRF_COLUMNS["ms2_mass_analyzer"]
        if ma_col in self.df.columns and row[ma_col]:
            params["ms2_mass_analyzer"] = self.parse_ontology_value(row[ma_col])

        # Scan window lower limit (DIA)
        sw_lower_col = SDRF_COLUMNS["scan_window_lower"]
        if sw_lower_col in self.df.columns and row[sw_lower_col]:
            params["scan_window_lower"] = self._parse_numeric_mz(row[sw_lower_col])

        # Scan window upper limit (DIA)
        sw_upper_col = SDRF_COLUMNS["scan_window_upper"]
        if sw_upper_col in self.df.columns and row[sw_upper_col]:
            params["scan_window_upper"] = self._parse_numeric_mz(row[sw_upper_col])

        # Isolation window width (DIA)
        iw_col = SDRF_COLUMNS["isolation_window_width"]
        if iw_col in self.df.columns and row[iw_col]:
            params["isolation_window_width"] = self._parse_numeric_mz(row[iw_col])

        return params

    @staticmethod
    def _parse_numeric_mz(value: str) -> Optional[float]:
        """Parse a numeric m/z value, stripping optional unit suffix like 'm/z'."""
        if not value or pd.isna(value):
            return None
        value = str(value).strip()
        # Accept "350", "350.5", "350 m/z", "350m/z"
        match = re.match(r"([\d.]+)\s*(?:m/z)?$", value, re.IGNORECASE)
        if match:
            try:
                return float(match.group(1))
            except ValueError:
                pass
        return None

    @staticmethod
    def parse_ontology_value(value: str) -> Dict[str, Optional[str]]:
        """
        Parse ontology-formatted value.

        Format: NT={term name};AC={accession}
        Example: NT=Q Exactive HF;AC=MS:1002523

        Args:
            value: The raw ontology string

        Returns:
            Dictionary with 'name', 'accession', and 'raw' keys
        """
        result = {"name": None, "accession": None, "raw": value}

        if not value or pd.isna(value):
            return result

        value = str(value).strip()

        # Parse NT=...
        nt_match = re.search(r"NT=([^;]+)", value, re.IGNORECASE)
        if nt_match:
            result["name"] = nt_match.group(1).strip()

        # Parse AC=...
        ac_match = re.search(r"AC=([^;]+)", value, re.IGNORECASE)
        if ac_match:
            result["accession"] = ac_match.group(1).strip()

        # If no NT/AC format, use the raw value as name
        if result["name"] is None and result["accession"] is None:
            result["name"] = value

        return result

    @staticmethod
    def parse_tolerance_value(value: str) -> Dict[str, Any]:
        """
        Parse tolerance value.

        Format: '10 ppm' or '0.02 Da'

        Args:
            value: The raw tolerance string

        Returns:
            Dictionary with 'value', 'unit', and 'raw' keys
        """
        result = {"value": None, "unit": None, "raw": value}

        if not value or pd.isna(value):
            return result

        value = str(value).strip()

        # Match number followed by unit
        match = re.match(r"([\d.]+)\s*(\w+)", value)
        if match:
            result["value"] = float(match.group(1))
            unit = match.group(2).lower()
            # Normalize unit names
            if unit == "da":
                result["unit"] = "Da"
            elif unit == "ppm":
                result["unit"] = "ppm"
            else:
                result["unit"] = unit

        return result
