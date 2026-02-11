"""
MassiveDownloader - Download files from MassIVE repository via FTP.
"""

import ftplib
import logging
import time
from pathlib import Path
from typing import Dict, Optional

from sdrf_refiner.downloader_base import BaseDownloader

logger = logging.getLogger(__name__)

MASSIVE_FTP_HOST = "massive-ftp.ucsd.edu"
MASSIVE_FTP_BASE = "/v03"

MAX_RETRIES = 3
RETRY_BACKOFF = 2  # seconds, doubled each retry


class MassiveDownloader(BaseDownloader):
    """
    Downloads raw files from MassIVE archive via FTP.

    On first use, builds an in-memory index of all files in the dataset
    by recursively walking the FTP directory tree.  Subsequent lookups
    are O(1) against this cached index.
    """

    def __init__(self, accession: str, download_dir: Path):
        super().__init__(accession, download_dir)
        self._file_index: Optional[Dict[str, str]] = None

    def _connect(self) -> ftplib.FTP:
        """Establish anonymous FTP connection with retry."""
        for attempt in range(MAX_RETRIES):
            try:
                ftp = ftplib.FTP(MASSIVE_FTP_HOST, timeout=60)
                ftp.login()  # anonymous
                logger.debug(f"Connected to {MASSIVE_FTP_HOST}")
                return ftp
            except ftplib.all_errors as e:
                if attempt < MAX_RETRIES - 1:
                    wait = RETRY_BACKOFF * (2 ** attempt)
                    logger.warning(
                        f"FTP connection failed (attempt {attempt + 1}): {e}. "
                        f"Retrying in {wait}s..."
                    )
                    time.sleep(wait)
                else:
                    raise ConnectionError(
                        f"Failed to connect to {MASSIVE_FTP_HOST} "
                        f"after {MAX_RETRIES} attempts: {e}"
                    ) from e

    def _build_file_index(self) -> Dict[str, str]:
        """
        Recursively list the FTP tree and build a filename -> path mapping.

        Returns:
            Dict mapping basename -> full FTP path
        """
        ftp = self._connect()
        base_path = f"{MASSIVE_FTP_BASE}/{self.accession}"
        index: Dict[str, str] = {}

        def _walk_mlsd(path: str) -> None:
            """Walk using MLSD (preferred, returns structured entries)."""
            try:
                for name, facts in ftp.mlsd(path):
                    if name in (".", ".."):
                        continue
                    full_path = f"{path}/{name}"
                    entry_type = facts.get("type", "").lower()
                    if entry_type == "dir":
                        _walk_mlsd(full_path)
                    elif entry_type == "file":
                        if name not in index:
                            index[name] = full_path
                        else:
                            logger.debug(
                                f"Duplicate filename '{name}' at {full_path} "
                                f"(keeping {index[name]})"
                            )
            except ftplib.error_perm:
                # MLSD not supported or permission denied, try NLST fallback
                _walk_nlst(path)

        def _walk_nlst(path: str) -> None:
            """Fallback walk using NLST (less structured)."""
            try:
                entries = ftp.nlst(path)
            except ftplib.all_errors as e:
                logger.warning(f"Could not list {path}: {e}")
                return

            for full_path in entries:
                name = full_path.rsplit("/", 1)[-1]
                if name in (".", ".."):
                    continue
                # Probe whether it's a directory by trying to CWD into it
                try:
                    ftp.cwd(full_path)
                    ftp.cwd("/")
                    _walk_nlst(full_path)
                except ftplib.error_perm:
                    # Not a directory — treat as file
                    if name not in index:
                        index[name] = full_path

        logger.info(
            f"Building file index for {self.accession} "
            f"(this may take a moment for large datasets)..."
        )
        _walk_mlsd(base_path)
        logger.info(f"Indexed {len(index)} files for {self.accession}")

        try:
            ftp.quit()
        except Exception:
            pass

        return index

    @property
    def file_index(self) -> Dict[str, str]:
        """Lazily build and cache the file index."""
        if self._file_index is None:
            self._file_index = self._build_file_index()
        return self._file_index

    def download_file_by_name(self, filename: str) -> Optional[Path]:
        """
        Download a single file from MassIVE by filename.

        Looks up the filename in the FTP index, then downloads via FTP RETR.

        Args:
            filename: The filename to download

        Returns:
            Path to downloaded file or None if failed
        """
        output_path = self.download_dir / filename

        # Check if already downloaded
        if output_path.exists():
            logger.info(f"File already exists: {filename}")
            return output_path

        # Look up in index
        ftp_path = self.file_index.get(filename)
        if ftp_path is None:
            logger.error(
                f"File '{filename}' not found in MassIVE dataset {self.accession}"
            )
            return None

        logger.info(f"Downloading {filename} from {self.accession} via FTP...")

        for attempt in range(MAX_RETRIES):
            try:
                ftp = self._connect()
                with open(output_path, "wb") as f:
                    ftp.retrbinary(
                        f"RETR {ftp_path}", f.write, blocksize=1024 * 1024
                    )
                try:
                    ftp.quit()
                except Exception:
                    pass

                if output_path.exists() and output_path.stat().st_size > 0:
                    logger.info(f"Downloaded: {filename}")
                    return output_path
                else:
                    logger.warning(f"Downloaded file is empty: {filename}")
                    if output_path.exists():
                        output_path.unlink()

            except ftplib.all_errors as e:
                if output_path.exists():
                    output_path.unlink()
                if attempt < MAX_RETRIES - 1:
                    wait = RETRY_BACKOFF * (2 ** attempt)
                    logger.warning(
                        f"FTP download failed (attempt {attempt + 1}): {e}. "
                        f"Retrying in {wait}s..."
                    )
                    time.sleep(wait)
                else:
                    logger.error(
                        f"Failed to download {filename} "
                        f"after {MAX_RETRIES} attempts: {e}"
                    )

        return None
