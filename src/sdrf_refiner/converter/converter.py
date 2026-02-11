"""
RawFileConverter - Convert vendor raw files to mzML format.
"""

import logging
import subprocess
from enum import Enum
from pathlib import Path
from typing import List, Optional

logger = logging.getLogger(__name__)


class VendorType(Enum):
    """Vendor types for raw files."""

    THERMO = "thermo"
    BRUKER = "bruker"
    SCIEX = "sciex"
    WATERS = "waters"
    UNKNOWN = "unknown"


class RawFileConverter:
    """
    Converts vendor-specific raw files to mzML format using appropriate tools.
    """

    # Extension to vendor mapping
    EXTENSION_VENDOR_MAP = {
        ".raw": VendorType.THERMO,
        ".d": VendorType.BRUKER,
        ".wiff": VendorType.SCIEX,
        ".wiff2": VendorType.SCIEX,
    }

    def __init__(self, mzml_dir: Path):
        """
        Initialize raw file converter.

        Args:
            mzml_dir: Directory to store converted mzML files
        """
        self.mzml_dir = Path(mzml_dir)
        self.mzml_dir.mkdir(parents=True, exist_ok=True)

    def detect_vendor(self, raw_file: Path) -> VendorType:
        """
        Detect vendor type from file extension.

        Args:
            raw_file: Path to raw file

        Returns:
            VendorType enum value
        """
        ext = raw_file.suffix.lower()
        return self.EXTENSION_VENDOR_MAP.get(ext, VendorType.UNKNOWN)

    def convert(self, raw_file: Path) -> Optional[Path]:
        """
        Convert a raw file to mzML format.

        Args:
            raw_file: Path to the raw file

        Returns:
            Path to the converted mzML file or None if conversion failed
        """
        raw_file = Path(raw_file)
        vendor = self.detect_vendor(raw_file)

        logger.info(f"Converting {raw_file.name} (vendor: {vendor.value})")

        try:
            if vendor == VendorType.THERMO:
                return self._convert_thermo(raw_file)
            elif vendor in [VendorType.BRUKER, VendorType.SCIEX, VendorType.WATERS]:
                return self._convert_msconvert(raw_file)
            else:
                logger.error(f"Unsupported file format: {raw_file.suffix}")
                return None
        except Exception as e:
            logger.error(f"Conversion failed for {raw_file.name}: {e}")
            return None

    def _convert_thermo(self, raw_file: Path) -> Path:
        """
        Convert Thermo RAW file using ThermoRawFileParser.

        Args:
            raw_file: Path to Thermo .raw file

        Returns:
            Path to converted mzML file
        """
        output_file = self.mzml_dir / f"{raw_file.stem}.mzML"

        if output_file.exists():
            logger.info(f"mzML already exists: {output_file.name}")
            return output_file

        # Try different command variations
        commands_to_try = [
            ["ThermoRawFileParser", f"-i={raw_file}", f"-o={self.mzml_dir}", "-f=2"],
            [
                "ThermoRawFileParser.sh",
                f"-i={raw_file}",
                f"-o={self.mzml_dir}",
                "-f=2",
            ],
            ["mono", "ThermoRawFileParser.exe", f"-i={raw_file}", f"-o={self.mzml_dir}", "-f=2"],
        ]

        last_error = None
        for cmd in commands_to_try:
            try:
                result = subprocess.run(
                    cmd, capture_output=True, text=True, check=True, timeout=600
                )
                logger.debug(result.stdout)

                if output_file.exists():
                    logger.info(f"Converted: {output_file.name}")
                    return output_file
            except FileNotFoundError:
                continue
            except subprocess.CalledProcessError as e:
                last_error = e
                logger.debug(f"Command failed: {e.stderr}")
                continue
            except subprocess.TimeoutExpired:
                last_error = RuntimeError("Conversion timed out after 10 minutes")
                continue

        raise RuntimeError(
            f"ThermoRawFileParser not found or failed. "
            f"Install via conda: conda install -c bioconda thermorawfileparser. "
            f"Last error: {last_error}"
        )

    def _convert_msconvert(self, raw_file: Path) -> Path:
        """
        Convert using ProteoWizard's msconvert.

        Args:
            raw_file: Path to raw file (Bruker, SCIEX, or Waters)

        Returns:
            Path to converted mzML file
        """
        output_file = self.mzml_dir / f"{raw_file.stem}.mzML"

        if output_file.exists():
            logger.info(f"mzML already exists: {output_file.name}")
            return output_file

        cmd = [
            "msconvert",
            str(raw_file),
            "-o",
            str(self.mzml_dir),
            "--mzML",
            "--64",
            "--filter",
            "peakPicking vendor msLevel=1-2",
        ]

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, check=True, timeout=600
            )
            logger.debug(result.stdout)

            if output_file.exists():
                logger.info(f"Converted: {output_file.name}")
                return output_file
            else:
                raise RuntimeError("mzML file not created after conversion")

        except FileNotFoundError:
            raise RuntimeError(
                "msconvert not found. Install ProteoWizard from: "
                "https://proteowizard.sourceforge.io/"
            )
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"msconvert failed: {e.stderr}")
        except subprocess.TimeoutExpired:
            raise RuntimeError("Conversion timed out after 10 minutes")

    def convert_all(self, raw_files: List[Path]) -> List[Path]:
        """
        Convert multiple raw files to mzML.

        Args:
            raw_files: List of paths to raw files

        Returns:
            List of paths to successfully converted mzML files
        """
        mzml_files = []

        for raw_file in raw_files:
            mzml = self.convert(raw_file)
            if mzml:
                mzml_files.append(mzml)

        logger.info(f"Converted {len(mzml_files)}/{len(raw_files)} files")
        return mzml_files

    def get_mzml_path(self, raw_file: Path) -> Path:
        """
        Get the expected mzML path for a raw file.

        Args:
            raw_file: Path to raw file

        Returns:
            Expected path for the mzML output
        """
        return self.mzml_dir / f"{raw_file.stem}.mzML"
