"""
TechSDRF - Validate and refine SDRF files using empirical MS data analysis.

This tool downloads raw MS files from PRIDE or MassIVE, analyzes them to detect
instrument parameters, and proposes refinements to SDRF metadata files.
"""

__version__ = "1.0.0"
__author__ = "PRIDE Team"

from sdrf_refiner.core.refiner import SDRFRefiner

__all__ = ["SDRFRefiner", "__version__"]
