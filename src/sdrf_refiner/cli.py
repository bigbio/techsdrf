#!/usr/bin/env python3
"""
SDRF Refiner CLI - Command-line interface for refining SDRF files
using empirical MS data analysis.
"""

import logging
import sys
from pathlib import Path
from typing import Optional, Tuple

import click

from sdrf_refiner import __version__
from sdrf_refiner.config import DEFAULT_NUM_FILES, SUPPORTED_FILE_TYPES
from sdrf_refiner.utils.logging import setup_logging


@click.group()
@click.version_option(version=__version__, prog_name="techsdrf")
def cli():
    """TechSDRF - Validate and refine SDRF files using MS data analysis."""
    pass


@cli.command("refine")
@click.option(
    "--accession",
    "-p",
    default=None,
    help="Dataset accession (e.g., PXD012345 for PRIDE, MSV000085836 for MassIVE)",
)
@click.option(
    "--pride-accession",
    hidden=True,
    default=None,
    help="[Deprecated] Use --accession/-p instead.",
)
@click.option(
    "--sdrf-file",
    "-s",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to existing SDRF file",
)
@click.option(
    "--num-files",
    "-n",
    default=DEFAULT_NUM_FILES,
    type=int,
    help=f"Number of raw files to download and analyze (default: {DEFAULT_NUM_FILES})",
)
@click.option(
    "--all-files",
    "-a",
    is_flag=True,
    help="Analyze all raw files (overrides --num-files)",
)
@click.option(
    "--output-sdrf",
    "-o",
    type=click.Path(path_type=Path),
    help="Output path for refined SDRF (default: <input>_refined.sdrf.tsv)",
)
@click.option(
    "--report-file",
    "-r",
    type=click.Path(path_type=Path),
    help="Output path for refinement report (default: <input>_report.txt)",
)
@click.option(
    "--work-dir",
    "-w",
    type=click.Path(path_type=Path),
    default=None,
    help="Working directory for downloads and processing (default: current directory)",
)
@click.option(
    "--keep-files",
    is_flag=True,
    help="Keep downloaded and converted files after processing. "
         "By default files are deleted after each file is analyzed to save disk space. "
         "Use this flag to retain files for resuming or debugging.",
)
@click.option(
    "--verbose",
    "-v",
    count=True,
    help="Increase verbosity (-v for INFO, -vv for DEBUG)",
)
@click.option(
    "--file-types",
    "-t",
    multiple=True,
    type=click.Choice(SUPPORTED_FILE_TYPES, case_sensitive=False),
    help="Restrict to specific file types (e.g., -t raw -t d)",
)
@click.option(
    "--tolerance-unit",
    "-u",
    type=click.Choice(["ppm", "Da"], case_sensitive=False),
    default=None,
    help="Preferred unit for tolerance values in refined SDRF (ppm or Da). "
         "If not set, uses the unit from the estimation method.",
)
@click.option(
    "--skip-ptm",
    is_flag=True,
    help="Skip PTM detection from spectra. "
         "By default, PTMs are detected and reported (but not written to SDRF).",
)
def refine(
    accession: Optional[str],
    pride_accession: Optional[str],
    sdrf_file: Path,
    num_files: int,
    all_files: bool,
    output_sdrf: Optional[Path],
    report_file: Optional[Path],
    work_dir: Optional[Path],
    keep_files: bool,
    verbose: int,
    file_types: Tuple[str, ...],
    tolerance_unit: Optional[str],
    skip_ptm: bool,
):
    """
    Refine an SDRF file by analyzing raw MS data from PRIDE or MassIVE.

    Downloads raw files, converts to mzML, analyzes instrument parameters,
    and proposes refinements to the SDRF based on detected values.

    \b
    Examples:
      techsdrf refine -p PXD012345 -s input.sdrf.tsv
      techsdrf refine -p MSV000085836 -s input.sdrf.tsv
      techsdrf refine --accession PXD012345 -s input.sdrf.tsv -n 5 -v
      techsdrf refine -p PXD012345 -s input.sdrf.tsv --all-files
      techsdrf refine -p PXD012345 -s input.sdrf.tsv -t raw --keep-files
    """
    # Set up logging
    setup_logging(verbose)
    logger = logging.getLogger("sdrf_refiner")

    # Handle deprecated --pride-accession alias
    effective_accession = accession or pride_accession
    if not effective_accession:
        click.secho("Error: --accession/-p is required.", fg="red")
        sys.exit(1)
    if pride_accession and not accession:
        logger.warning("--pride-accession is deprecated, use --accession/-p instead")

    # Validate accession format early
    from sdrf_refiner.downloader_factory import detect_repository

    try:
        repo = detect_repository(effective_accession)
        click.echo(f"Repository: {repo} ({effective_accession.upper()})")
    except ValueError as e:
        click.secho(str(e), fg="red")
        sys.exit(1)

    try:
        from sdrf_refiner.core.refiner import SDRFRefiner

        # Use None for num_files when --all-files is set
        effective_num_files = None if all_files else num_files

        refiner = SDRFRefiner(
            accession=effective_accession,
            sdrf_file=sdrf_file,
            num_files=effective_num_files,
            output_sdrf=output_sdrf,
            report_file=report_file,
            work_dir=work_dir,
            keep_files=keep_files,
            file_types=list(file_types) if file_types else None,
            tolerance_unit=tolerance_unit,
            skip_ptm=skip_ptm,
        )

        result = refiner.run()

        if result.success:
            click.echo()
            click.secho(f"Refinement complete!", fg="green", bold=True)
            click.echo(f"  Refined SDRF: {result.output_sdrf}")
            click.echo(f"  Report: {result.report_file}")
            click.echo(f"  Refinements applied: {result.num_refinements}")
            sys.exit(0)
        else:
            click.secho(f"Refinement failed: {result.error_message}", fg="red")
            sys.exit(1)

    except Exception as e:
        click.secho(f"Error: {str(e)}", fg="red")
        if verbose >= 2:
            import traceback

            traceback.print_exc()
        sys.exit(1)


@cli.command("validate")
@click.option(
    "--sdrf-file",
    "-s",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to SDRF file to validate",
)
@click.option(
    "--verbose",
    "-v",
    count=True,
    help="Increase verbosity",
)
def validate(sdrf_file: Path, verbose: int):
    """
    Validate an SDRF file using sdrf-pipelines.

    \b
    Example:
      techsdrf validate -s input.sdrf.tsv
    """
    setup_logging(verbose)

    from sdrf_refiner.sdrf.reader import SDRFReader

    reader = SDRFReader(sdrf_file)
    errors = reader.validate()

    if not errors:
        click.secho("SDRF validation passed!", fg="green")
        sys.exit(0)
    else:
        click.secho(f"Validation found {len(errors)} issues:", fg="yellow")
        for error in errors:
            click.echo(f"  - {error}")
        sys.exit(1)


@cli.command("info")
@click.option(
    "--sdrf-file",
    "-s",
    required=True,
    type=click.Path(exists=True, path_type=Path),
    help="Path to SDRF file",
)
def info(sdrf_file: Path):
    """
    Display information about an SDRF file.

    Shows the parameters currently defined in the SDRF.

    \b
    Example:
      techsdrf info -s input.sdrf.tsv
    """
    from sdrf_refiner.sdrf.reader import SDRFReader

    reader = SDRFReader(sdrf_file)
    reader.read()

    click.echo(f"SDRF File: {sdrf_file}")
    click.echo(f"Rows: {len(reader.df)}")
    click.echo(f"Columns: {len(reader.df.columns)}")
    click.echo()

    params = reader.get_all_parameters()

    click.echo("Current Parameters:")
    click.echo("-" * 40)

    if not params:
        click.echo("  (no instrument parameters found)")
    else:
        for key, value in params.items():
            if value:
                if isinstance(value, dict):
                    name = value.get("name") or value.get("raw", "")
                    val = value.get("value")
                    unit = value.get("unit", "")
                    if val is not None:
                        click.echo(f"  {key}: {val} {unit}")
                    elif name:
                        click.echo(f"  {key}: {name}")
                else:
                    click.echo(f"  {key}: {value}")
            else:
                click.echo(f"  {key}: (not specified)")

    click.echo()

    data_files = reader.get_data_files()
    if data_files:
        click.echo(f"Data files referenced: {len(data_files)}")
        for f in data_files[:5]:
            click.echo(f"  - {f}")
        if len(data_files) > 5:
            click.echo(f"  ... and {len(data_files) - 5} more")


def main():
    """Entry point for the CLI."""
    cli()


if __name__ == "__main__":
    main()
