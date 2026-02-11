"""Working directory management for downloads and processing."""

import logging
import shutil
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


class TempManager:
    """
    Manages working directory for downloads and converted files.

    Uses the provided work_dir or the current working directory (cwd).
    Never uses system temp; files stay in the chosen folder.
    """

    def __init__(self, work_dir: Optional[Path] = None, keep_files: bool = False):
        """
        Initialize working directory manager.

        Args:
            work_dir: Working directory for downloads and mzML. If None, uses cwd.
            keep_files: If True, keep downloaded/converted files after processing.
                        Useful for resuming after pipeline changes.
        """
        self.keep_files = keep_files
        self._is_temp = False  # Never use system temp

        if work_dir:
            self.work_dir = Path(work_dir).resolve()
        else:
            self.work_dir = Path.cwd()

        self.work_dir.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        self.download_dir = self.work_dir / "downloads"
        self.mzml_dir = self.work_dir / "mzml"
        self.download_dir.mkdir(exist_ok=True)
        self.mzml_dir.mkdir(exist_ok=True)

        logger.info(f"Working directory: {self.work_dir} (downloads and mzML stored here)")

    def remove_file(self, file_path: Path) -> None:
        """
        Remove a single file or directory.

        Args:
            file_path: Path to the file or directory to remove
        """
        if self.keep_files:
            return
        try:
            file_path = Path(file_path)
            if file_path.is_dir():
                shutil.rmtree(file_path)
                logger.debug(f"Removed directory: {file_path.name}")
            elif file_path.exists():
                file_path.unlink()
                logger.debug(f"Removed file: {file_path.name}")
        except OSError as e:
            logger.warning(f"Could not remove {file_path}: {e}")

    def cleanup(self):
        """Clean up only when using temp (not used anymore)."""
        if self.keep_files:
            logger.info(f"Keeping files at: {self.work_dir}")
            return
        # No cleanup for user-specified or cwd work_dir; files remain

    def __enter__(self):
        """Context manager entry."""
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit with cleanup."""
        self.cleanup()
        return False

    def __del__(self):
        """Destructor with cleanup."""
        try:
            self.cleanup()
        except Exception:
            pass  # Ignore errors during garbage collection
