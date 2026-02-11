"""
BaseDownloader - Abstract base class for MS file downloaders.

Provides shared logic for filename resolution (filtering, deduplication,
limiting) and defines the interface that repository-specific downloaders
must implement.
"""

import logging
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Optional

logger = logging.getLogger(__name__)


class BaseDownloader(ABC):
    """
    Abstract base for downloading MS data files from a repository.

    Subclasses must implement ``download_file_by_name``.  Shared helpers
    (``resolve_filenames``, ``download_files_from_sdrf``) are provided here.
    """

    # File extensions we can handle (raw vendor files + already-converted mzML)
    EXTENSION_MAP: Dict[str, str] = {
        ".raw": "thermo",
        ".d": "bruker",
        ".wiff": "sciex",
        ".wiff2": "sciex",
        ".mzml": "mzml",
    }

    def __init__(self, accession: str, download_dir: Path):
        self.accession = accession.upper()
        self.download_dir = Path(download_dir)
        self.download_dir.mkdir(parents=True, exist_ok=True)

    def resolve_filenames(
        self,
        filenames: List[str],
        num_files: Optional[int] = None,
        file_types: Optional[List[str]] = None,
    ) -> List[str]:
        """
        Resolve the list of filenames to process from SDRF data file entries.

        Filters by supported extensions, deduplicates, and limits to
        ``num_files``.  Does **not** download anything.

        Args:
            filenames: List of filenames from SDRF comment[data file] column
            num_files: Maximum number of files to return (None = all)
            file_types: Optional filter for file types (e.g., ['raw', 'd'])

        Returns:
            Ordered list of unique filenames to process
        """
        if file_types:
            extensions = [f".{ft.lower().lstrip('.')}" for ft in file_types]
        else:
            extensions = list(self.EXTENSION_MAP.keys())

        # Filter filenames by extension
        raw_filenames = []
        for filename in filenames:
            ext = Path(filename).suffix.lower()
            if ext in extensions:
                raw_filenames.append(filename)

        # Remove duplicates while preserving order
        seen: set = set()
        unique_filenames: List[str] = []
        for f in raw_filenames:
            if f not in seen:
                seen.add(f)
                unique_filenames.append(f)

        logger.info(f"Found {len(unique_filenames)} unique raw files in SDRF")

        # Limit number of files if specified
        if num_files is not None and len(unique_filenames) > num_files:
            unique_filenames = unique_filenames[:num_files]
            logger.info(f"Limited to {num_files} files")

        return unique_filenames

    @abstractmethod
    def download_file_by_name(self, filename: str) -> Optional[Path]:
        """
        Download a single file by name.

        Args:
            filename: The filename to download

        Returns:
            Path to downloaded file or None if failed
        """
        ...

    def download_files_from_sdrf(
        self,
        filenames: List[str],
        num_files: Optional[int] = None,
        file_types: Optional[List[str]] = None,
    ) -> List[Path]:
        """
        Download files based on filenames from SDRF.

        Args:
            filenames: List of filenames from SDRF comment[data file] column
            num_files: Number of files to download (None = all files)
            file_types: Optional filter for file types (e.g., ['raw', 'd'])

        Returns:
            List of paths to downloaded files
        """
        unique_filenames = self.resolve_filenames(filenames, num_files, file_types)

        downloaded = []
        for filename in unique_filenames:
            path = self.download_file_by_name(filename)
            if path:
                downloaded.append(path)
            else:
                logger.warning(f"Failed to download {filename}")

        logger.info(f"Successfully downloaded {len(downloaded)}/{len(unique_filenames)} files")
        return downloaded
