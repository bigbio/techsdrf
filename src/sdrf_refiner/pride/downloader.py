"""
PrideDownloader - Download files from PRIDE using pridepy or direct HTTP.
"""

import logging
from pathlib import Path
from typing import Any, Dict, List, Optional

import requests

from sdrf_refiner.downloader_base import BaseDownloader

logger = logging.getLogger(__name__)

PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v3"
FTP_BASE = "ftp://ftp.pride.ebi.ac.uk/"
GLOBUS_HTTP_BASE = "https://g-a8b222.dd271.03c0.data.globus.org/"


class PrideDownloader(BaseDownloader):
    """
    Downloads raw files from PRIDE archive using pridepy.
    """

    def __init__(self, pride_accession: str, download_dir: Path):
        """
        Initialize PRIDE downloader.

        Args:
            pride_accession: PRIDE accession (e.g., PXD012345)
            download_dir: Directory to download files to
        """
        super().__init__(pride_accession, download_dir)
        self.pride_accession = self.accession  # backward-compat alias

    # Protocol order: try direct HTTP (no timeout) first, then pridepy
    _DOWNLOAD_PROTOCOLS = ("ftp", "globus")

    def _fetch_file_metadata(self, filename: str) -> Optional[Dict[str, Any]]:
        """Fetch file metadata from PRIDE API (no timeout)."""
        url = f"{PRIDE_API_BASE}/projects/{self.pride_accession}/files/all"
        try:
            resp = requests.get(
                url, headers={"Accept": "application/json"}, stream=False, timeout=60
            )
            resp.raise_for_status()
            files = resp.json()
            for f in files:
                if f.get("fileName") == filename:
                    return f
        except Exception as e:
            logger.debug(f"PRIDE API request failed: {e}")
        return None

    def _download_via_http_no_timeout(self, http_url: str, output_path: Path) -> bool:
        """
        Download file via HTTP with no timeout (for large files).
        Uses streaming with no read/connect timeout.
        """
        try:
            session = requests.Session()
            session.headers["User-Agent"] = "techsdrf/1.0"
            # No timeout - allow arbitrarily long downloads for large RAW files
            with session.get(http_url, stream=True, timeout=None) as r:
                r.raise_for_status()
                total = int(r.headers.get("content-length", 0))
                with open(output_path, "wb") as f:
                    downloaded = 0
                    last_log_mb = 0
                    for chunk in r.iter_content(chunk_size=1024 * 1024):
                        if chunk:
                            f.write(chunk)
                            downloaded += len(chunk)
                            mb = downloaded // (1024 * 1024)
                            if mb >= last_log_mb + 50:
                                last_log_mb = mb
                                pct = 100 * downloaded / total if total else 0
                                logger.info(f"  Progress: {downloaded / (1024**3):.1f} GB ({pct:.0f}%)")
            return output_path.exists()
        except Exception as e:
            logger.warning(f"HTTP download failed: {e}")
            if output_path.exists():
                output_path.unlink()
            return False

    def download_file_by_name(self, filename: str) -> Optional[Path]:
        """
        Download a single file from PRIDE by filename.

        Uses direct HTTP (no timeout) first, then falls back to pridepy (FTP/globus).

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

        logger.info(f"Downloading {filename} from {self.pride_accession}...")

        # 1. Try direct HTTP (no timeout) via Globus mirror
        meta = self._fetch_file_metadata(filename)
        if meta and "publicFileLocations" in meta:
            for loc in meta["publicFileLocations"]:
                val = loc.get("value", "")
                if val.startswith(FTP_BASE):
                    http_url = val.replace(FTP_BASE, GLOBUS_HTTP_BASE)
                    logger.info(f"Using HTTP download (no timeout)")
                    if self._download_via_http_no_timeout(http_url, output_path):
                        logger.info(f"Downloaded: {filename} (via HTTP)")
                        return output_path
                    break

        # 2. Fall back to pridepy (FTP, globus) - may have timeouts
        try:
            from pridepy.files.files import Files

            pride_files = Files()
            last_error: Optional[Exception] = None

            for protocol in self._DOWNLOAD_PROTOCOLS:
                try:
                    pride_files.download_file_by_name(
                        accession=self.pride_accession,
                        file_name=filename,
                        output_folder=str(self.download_dir),
                        skip_if_downloaded_already=True,
                        protocol=protocol,
                        username=None,
                        password=None,
                        aspera_maximum_bandwidth="100M",
                        checksum_check=False,
                    )
                    if output_path.exists():
                        logger.info(f"Downloaded: {filename} (via {protocol})")
                        return output_path
                except Exception as e:
                    last_error = e
                    logger.warning(f"Download via {protocol} failed: {e}. Trying next protocol...")

            logger.error(f"Failed to download {filename} after trying all protocols: {last_error}")
            return None

        except ImportError:
            logger.error("pridepy not installed. Install with: pip install pridepy")
            return None

    # Keep backward compatibility with old API-based method
    def download_files(
        self, num_files: Optional[int], file_types: Optional[List[str]] = None
    ) -> List[Path]:
        """
        Legacy method - requires filenames to be provided via download_files_from_sdrf.

        This method is kept for compatibility but should not be used directly.
        Use download_files_from_sdrf instead.
        """
        logger.warning(
            "download_files() called without filenames. "
            "Use download_files_from_sdrf() with filenames from SDRF."
        )
        return []
