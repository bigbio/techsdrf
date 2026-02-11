"""
Factory for creating the appropriate downloader based on accession prefix.
"""

import re
from pathlib import Path

from sdrf_refiner.downloader_base import BaseDownloader

# Accession patterns
PXD_PATTERN = re.compile(r"^PXD\d+$", re.IGNORECASE)
MSV_PATTERN = re.compile(r"^MSV\d{9}$", re.IGNORECASE)


def create_downloader(accession: str, download_dir: Path) -> BaseDownloader:
    """
    Create the appropriate downloader based on accession format.

    Args:
        accession: Repository accession (PXD... for PRIDE, MSV... for MassIVE)
        download_dir: Directory to download files into

    Returns:
        A BaseDownloader subclass instance

    Raises:
        ValueError: If accession format is not recognized
    """
    accession_upper = accession.strip().upper()

    if PXD_PATTERN.match(accession_upper):
        from sdrf_refiner.pride.downloader import PrideDownloader

        return PrideDownloader(accession_upper, download_dir)

    elif MSV_PATTERN.match(accession_upper):
        from sdrf_refiner.massive.downloader import MassiveDownloader

        return MassiveDownloader(accession_upper, download_dir)

    else:
        raise ValueError(
            f"Unrecognized accession format: '{accession}'. "
            f"Expected PXD (PRIDE) or MSV (MassIVE) prefix."
        )


def detect_repository(accession: str) -> str:
    """
    Detect repository type from accession string.

    Returns:
        'PRIDE' or 'MassIVE'

    Raises:
        ValueError: If format not recognized
    """
    accession_upper = accession.strip().upper()
    if PXD_PATTERN.match(accession_upper):
        return "PRIDE"
    elif MSV_PATTERN.match(accession_upper):
        return "MassIVE"
    else:
        raise ValueError(
            f"Unrecognized accession format: '{accession}'. "
            f"Expected PXD (PRIDE) or MSV (MassIVE) prefix."
        )
