"""
SDRFWriter - Write refined SDRF files.
"""

import logging
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

from sdrf_refiner.config import SDRF_COLUMNS, format_ontology_value
from sdrf_refiner.sdrf.comparer import ParameterComparison, ComparisonStatus

logger = logging.getLogger(__name__)

# Regex matching the pandas duplicate-column suffix (.1, .2, …)
_PANDAS_DUP_SUFFIX = re.compile(r"\.\d+$")


class SDRFWriter:
    """
    Writes refined SDRF files with updated parameter values.

    Handles SDRF files with duplicate column names properly.
    """

    @staticmethod
    def _strip_pandas_suffixes(pandas_columns: List[str]) -> List[str]:
        """
        Remove the ``.1``, ``.2``, … suffixes that pandas appends to duplicate
        column names so the header written to disk matches the original file.

        Only strips a suffix when the base name already appeared earlier in the
        list (i.e. it really is a pandas-generated duplicate, not a column whose
        real name ends with a digit).
        """
        seen: Dict[str, int] = {}
        clean: List[str] = []
        for col in pandas_columns:
            m = _PANDAS_DUP_SUFFIX.search(col)
            if m:
                base = col[: m.start()]
                if base in seen:
                    clean.append(base)
                    seen[base] += 1
                    continue
            seen[col] = seen.get(col, 0) + 1
            clean.append(col)
        return clean

    def __init__(self, original_df: pd.DataFrame, original_columns: Optional[List[str]] = None):
        """
        Initialize SDRF writer.

        Args:
            original_df: The original SDRF DataFrame to modify
            original_columns: Original column names (preserving duplicates).
                If *None*, the clean names are derived from the DataFrame
                columns by stripping pandas duplicate suffixes (``.1``, …).
        """
        self.df = original_df.copy()
        # Positional mapping: pandas may rename duplicate columns with .1, .2
        # suffixes.  We keep the pandas column names so we can look up values
        # by position rather than by (ambiguous) original name.
        self._pandas_columns: List[str] = list(self.df.columns)
        # _original_columns holds the *clean* names (no .1/.2 suffixes) that
        # will be written as the TSV header.
        if original_columns is not None:
            self._original_columns = list(original_columns)
        else:
            self._original_columns = self._strip_pandas_suffixes(self._pandas_columns)
        self.changes: List[Dict[str, Any]] = []  # Track changes for reporting

    def apply_refinements(self, refinements: List[ParameterComparison]) -> int:
        """
        Apply refinements to the SDRF DataFrame.

        Args:
            refinements: List of parameter comparisons with recommendations

        Returns:
            Number of refinements applied
        """
        applied = 0
        self.changes = []

        for refinement in refinements:
            if refinement.status not in [
                ComparisonStatus.MISMATCH,
                ComparisonStatus.MISSING_SDRF,
                ComparisonStatus.IMPROVED,
            ]:
                continue

            # Never apply tolerance changes that are heuristic-only
            if refinement.details.get("heuristic_only"):
                logger.info(
                    f"Skipping heuristic-only tolerance for {refinement.parameter_name} "
                    f"(no empirical estimation available)"
                )
                continue

            col_name = self._param_to_column(refinement.parameter_name)
            if not col_name:
                logger.warning(f"Unknown parameter: {refinement.parameter_name}")
                continue

            new_value = self._format_value(refinement)
            if not new_value:
                continue

            # Get old value for reporting
            old_value = ""
            if col_name in self.df.columns:
                old_value = self.df[col_name].iloc[0] if len(self.df) > 0 else ""
            else:
                # Add column if missing
                self.df[col_name] = ""
                self._original_columns.append(col_name)
                self._pandas_columns.append(col_name)
                logger.info(f"Added new column: {col_name}")

            # Update all rows with the new value
            self.df[col_name] = new_value

            # Track the change for reporting
            self.changes.append({
                "column": col_name,
                "old_value": old_value if old_value else "(not specified)",
                "new_value": new_value,
            })

            logger.info(f"{col_name}: {old_value if old_value else '(not specified)'} -> {new_value}")
            applied += 1

        return applied

    def get_changes(self) -> List[Dict[str, Any]]:
        """
        Get list of changes made.

        Returns:
            List of dicts with 'column', 'old_value', 'new_value'
        """
        return self.changes

    def write(self, output_path: Path) -> None:
        """
        Write the refined SDRF to a file.

        Preserves original column names including duplicates.  Uses positional
        indexing so that duplicate column names (e.g. multiple
        ``comment[modification parameters]``) each write their own data instead
        of all resolving to the first column with that name.

        Args:
            output_path: Path for the output file
        """
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # Write manually to preserve duplicate column names
        with open(output_path, "w", encoding="utf-8") as f:
            # Write header with original column names
            f.write("\t".join(self._original_columns) + "\n")

            # Write data rows using positional indexing so that duplicate
            # column names each map to their own pandas column (.1, .2, …)
            for idx in range(len(self.df)):
                row_values = []
                for col_pos in range(len(self._original_columns)):
                    if col_pos < len(self._pandas_columns):
                        pandas_col = self._pandas_columns[col_pos]
                        val = self.df[pandas_col].iloc[idx]
                    else:
                        # Fallback for columns added after init
                        val = self._get_column_value(
                            self._original_columns[col_pos], idx
                        )
                    row_values.append(str(val) if val and not pd.isna(val) else "")
                f.write("\t".join(row_values) + "\n")

        logger.info(f"Written refined SDRF to: {output_path}")

    def _get_column_value(self, col_name: str, row_idx: int) -> str:
        """Get value for a column, handling pandas renamed duplicates."""
        # Direct match
        if col_name in self.df.columns:
            return self.df[col_name].iloc[row_idx]

        # Try to find pandas-renamed column (e.g., "col.1", "col.2")
        for df_col in self.df.columns:
            if df_col.startswith(col_name + "."):
                return self.df[df_col].iloc[row_idx]

        return ""

    @staticmethod
    def _param_to_column(param_name: str) -> Optional[str]:
        """Map parameter name to SDRF column name."""
        mapping = {
            "instrument": SDRF_COLUMNS["instrument"],
            "dissociation_method": SDRF_COLUMNS["dissociation"],
            "ms3_dissociation_method": SDRF_COLUMNS["ms3_dissociation"],
            "precursor_mass_tolerance": SDRF_COLUMNS["precursor_tolerance"],
            "fragment_mass_tolerance": SDRF_COLUMNS["fragment_tolerance"],
            "ms_min_charge": SDRF_COLUMNS["ms_min_charge"],
            "ms_max_charge": SDRF_COLUMNS["ms_max_charge"],
            "collision_energy": SDRF_COLUMNS["collision_energy"],
            "ms2_mass_analyzer": SDRF_COLUMNS["ms2_mass_analyzer"],
            "scan_window_lower": SDRF_COLUMNS["scan_window_lower"],
            "scan_window_upper": SDRF_COLUMNS["scan_window_upper"],
            "isolation_window_width": SDRF_COLUMNS["isolation_window_width"],
        }
        return mapping.get(param_name)

    # Parameters whose detected values contain ontology dicts (name + accession)
    _ONTOLOGY_PARAMS = frozenset({
        "instrument",
        "dissociation_method",
        "ms3_dissociation_method",
        "ms2_mass_analyzer",
    })

    @staticmethod
    def _format_value(refinement: ParameterComparison) -> Optional[str]:
        """
        Format detected value for SDRF output.

        Ontology-valued columns use ``AC=<accession>;NT=<name>`` format.
        Tolerance columns use ``<value> <unit>``.
        Charge columns are plain integers.
        Collision energy is a plain string.

        Args:
            refinement: The parameter comparison with detected value

        Returns:
            Formatted string for SDRF or None
        """
        detected = refinement.detected_value
        if not detected:
            return None

        param = refinement.parameter_name

        # Ontology-valued columns → AC=…;NT=…
        if param in SDRFWriter._ONTOLOGY_PARAMS:
            if isinstance(detected, dict):
                name = detected.get("name", "")
                accession = detected.get("accession", "")
                return format_ontology_value(name, accession) or None
            return None

        if param in ["ms_min_charge", "ms_max_charge"]:
            # Charge values are simple integers
            if isinstance(detected, int):
                return str(detected)
            return None

        if "tolerance" in param:
            if isinstance(detected, dict):
                value = detected.get("value")
                unit = detected.get("unit", "ppm")
                if value is not None:
                    # Round appropriately based on unit
                    if unit == "ppm":
                        return f"{int(round(value))} ppm"
                    else:
                        return f"{value:.4f} {unit}".rstrip("0").rstrip(".")
            return None

        if param == "collision_energy":
            if isinstance(detected, str) and detected.strip():
                return detected.strip()
            return None

        # Numeric m/z values (scan windows, isolation width) → plain number
        if param in ("scan_window_lower", "scan_window_upper", "isolation_window_width"):
            if isinstance(detected, (int, float)):
                # Use integer format when the value is whole, otherwise 2 decimals
                return f"{detected:g}"
            return None

        return None
