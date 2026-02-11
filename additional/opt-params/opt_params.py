import logging
import math
import random
import re
import uuid
from collections import Counter
from typing import Tuple, List, Dict
import json
import os
from sagepy.core import EnzymeBuilder, SageSearchConfiguration, validate_mods, validate_var_mods, \
    Scorer, RawSpectrum, SpectrumProcessor, Precursor, Tolerance, IonType
from sagepy.qfdr.tdc import target_decoy_competition_pandas
from sagepy.utility import mean_ppm, median_ppm
from sagepy.core.ml.pep import calculate_pep
from sagepy.rescore.lda import rescore_lda

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pandas import DataFrame
import click
from pyteomics import mzml, mass
import concurrent.futures as cf

logging.basicConfig(level=logging.INFO)

folder_plots_uui = ""
folder_sperator = os.path.sep

CONTEXT_SETTINGS = dict(help_option_names=["-h", "--help"])

indexed_db = None


class SageRunner:
    def __init__(self):
        self.queries = []

    def generate_indexed_db(self, fasta_path, param_data):
        # if sage_config_file is None:
        #     raise ValueError("The sage config file is required.")
        #
        # with open(sage_config_file) as f:
        #     param_data = json.load(f)
        cleave_at = param_data["database"]["enzyme"]["cleave_at"]
        missed_cleavages = param_data["database"]["enzyme"]["missed_cleavages"]
        min_len = param_data["database"]["enzyme"]["min_len"]
        max_len = param_data["database"]["enzyme"]["max_len"]
        restrict = param_data["database"]["enzyme"]["restrict"]
        static_mods = param_data["database"]["static_mods"]
        variable_mods = param_data["database"]["variable_mods"]
        fragment_min_mz = param_data["database"]["fragment_min_mz"]
        fragment_max_mz = param_data["database"]["fragment_max_mz"]
        peptide_min_mass = param_data["database"]["peptide_min_mass"]
        peptide_max_mass = param_data["database"]["peptide_max_mass"]
        ion_kinds = param_data["database"]["ion_kinds"]
        min_ion_index = param_data["database"]["min_ion_index"]
        max_variable_mods = param_data["database"]["max_variable_mods"]
        generate_decoys = param_data["database"]["generate_decoys"]
        decoy_tag = param_data["database"]["decoy_tag"]

        # configure a trypsin-like digestor of fasta files
        enzyme_builder = EnzymeBuilder(
            missed_cleavages=missed_cleavages,
            min_len=min_len,
            max_len=max_len,
            cleave_at=cleave_at,
            restrict=restrict,
            c_terminal=True
        )

        with open(fasta_path, 'r') as infile:
            fasta = infile.read()

        # generate IonType class
        ion_kinds = [IonType(ion).get_py_ptr() for ion in ion_kinds]

        # set-up a config for a sage-database
        sage_config = SageSearchConfiguration(
            fasta=fasta,
            static_mods=static_mods,
            variable_mods=variable_mods,
            enzyme_builder=enzyme_builder,
            fragment_min_mz=fragment_min_mz,
            fragment_max_mz=fragment_max_mz,
            peptide_min_mass=peptide_min_mass,
            peptide_max_mass=peptide_max_mass,
            ion_kinds=ion_kinds,
            generate_decoys=generate_decoys,
            decoy_tag=decoy_tag,
            min_ion_index=min_ion_index,
            max_variable_mods=max_variable_mods,
            bucket_size=int(np.power(2, 14))
        )

        # generate the database for searching against
        indexed_db = sage_config.generate_indexed_database()
        return indexed_db

    def peptide_spectrum_match_list_to_pandas(self,
                                              psms, indexed_db,
                                              re_score: bool = False,
                                              use_sequence_as_match_idx: bool = True) -> pd.DataFrame:
        """Convert a list of peptide spectrum matches to a pandas dataframe

        Args:
            psms (List[PeptideSpectrumMatch]): The peptide spectrum matches
            re_score (bool, optional): Should re-score be used. Defaults to False.
            use_sequence_as_match_idx (bool, optional): Should the sequence be used as the match index. Defaults to True.

        Returns:
            pd.DataFrame: The pandas dataframe
        """
        row_list = []

        for match in psms:
            if match.retention_time_predicted is not None and match.projected_rt is not None:
                delta_rt = match.retention_time_predicted - match.projected_rt
            else:
                delta_rt = None

            if re_score:
                score = match.re_score
            else:
                score = match.hyper_score

            if use_sequence_as_match_idx:
                match_idx = match.sequence
            else:
                match_idx = str(match.peptide_idx)

            if match.inverse_mobility_predicted is not None:
                delta_ims = match.inverse_mobility_predicted - match.inverse_mobility_observed
            else:
                delta_ims = None

            if match.beta_score is not None:
                beta_score = match.beta_score
            else:
                beta_score = None

            peptide = indexed_db[match.peptide_idx]
            position = peptide.position
            modifications = peptide.modifications
            precursor_ppm = abs(match.mono_mass_observed - match.mono_mass_calculated - match.isotope_error) * 2e6 / (
                    match.mono_mass_observed + match.mono_mass_calculated - match.isotope_error)

            row_list.append({
                "spec_idx": match.spec_idx,
                "match_idx": match_idx,
                "proteins": ";".join(match.proteins),
                "decoy": match.decoy,
                "score": score,
                "re_score": match.re_score,
                "hyperscore": match.hyper_score,
                "rank": match.rank,
                "mono_mz_calculated": match.mono_mz_calculated,
                "mono_mass_observed": match.mono_mass_observed,
                "mono_mass_calculated": match.mono_mass_calculated,
                "delta_mass": match.mono_mass_calculated - match.mono_mass_observed,
                "isotope_error": match.isotope_error,
                "fragment_ppm": match.average_ppm,
                "precursor_ppm": precursor_ppm,
                "delta_next": match.delta_next,
                "delta_best": match.delta_best,
                "matched_peaks": match.matched_peaks,
                "longest_b": match.longest_b,
                "longest_y": match.longest_y,
                "longest_y_pct": match.longest_y_pct,
                "missed_cleavages": match.missed_cleavages,
                "matched_intensity_pct": match.matched_intensity_pct,
                "scored_candidates": match.scored_candidates,
                "poisson": match.poisson,
                "modifications": modifications,
                "charge": match.charge,
                "retention_time_observed": match.retention_time_observed,
                "retention_time_predicted": match.retention_time_predicted,
                "delta_rt": delta_rt,
                "inverse_mobility_observed": match.inverse_mobility_observed,
                "inverse_mobility_predicted": match.inverse_mobility_predicted,
                "delta_ims": delta_ims,
                "intensity_ms1": match.intensity_ms1,
                "intensity_ms2": match.intensity_ms2,
                "q_value": match.q_value,
                "collision_energy": match.collision_energy,
                "cosine_similarity": match.cosine_similarity,
                "projected_rt": match.projected_rt,
                "beta_score": beta_score
            })

        return pd.DataFrame(row_list)

    def parse_mzml(self, spec_file, file_id, fragment_min_mz, fragment_max_mz, deisotope_bool):
        raw_spectrums = mzml.MzML(spec_file)

        for spectrum in raw_spectrums:
            if spectrum["ms level"] == 1:
                continue
            precursor = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]
            precursor_mz = precursor["selected ion m/z"]
            precursor_charge = precursor["charge state"]
            precursor_intensity = precursor.get('peak intensity', None)
            collision_energy = spectrum['precursorList']['precursor'][0]['activation']['collision energy']
            injection_time = spectrum['scanList']['scan'][0]['ion injection time']
            retention_time = spectrum['scanList']['scan'][0]['scan start time']

            # Extract the relevant metadata
            isolation_window = spectrum['precursorList']['precursor'][0]['isolationWindow']

            lower = isolation_window['isolation window target m/z'] - isolation_window[
                'isolation window lower offset']
            upper = isolation_window['isolation window target m/z'] + isolation_window[
                'isolation window upper offset']

            mz = spectrum["m/z array"].astype(np.float32)
            intensity = spectrum["intensity array"].astype(np.float32)
            # set selection window bounds
            tolerance = Tolerance(da=(lower, upper))
            precursor = Precursor(
                charge=precursor_charge,
                mz=precursor_mz,
                intensity=precursor_intensity,
                isolation_window=tolerance,
                collision_energy=collision_energy
            )
            scan = re.findall(r"scan=(\d+)", spectrum["id"])[0]
            raw_spectrum = RawSpectrum(
                file_id=file_id,
                spec_id=os.path.splitext(os.path.basename(spec_file))[0] + "." + str(scan) + "." + str(
                    scan) + "." + str(precursor_charge),
                total_ion_current=spectrum["total ion current"],
                ion_injection_time=injection_time,
                scan_start_time=retention_time,
                precursors=[precursor],
                mz=mz,
                intensity=intensity
            )

            spec_processor = SpectrumProcessor(take_top_n=150, min_fragment_mz=fragment_min_mz,
                                               max_fragment_mz=fragment_max_mz,
                                               deisotope=deisotope_bool)

            query = spec_processor.process(raw_spectrum)
            self.queries.append(query)

        # return querys

    def run_sage(self, fragment_tolerance: int = 0, precursor_tolerance: int = 0, fragment_type: str = "ppm",
                 mzml_files: list = [], sage_config_file: str = None, fasta_file: str = None,
                 use_file_values: bool = True, skip_spectra_preprocess: bool = False) -> DataFrame:
        if sage_config_file is None:
            raise ValueError("The sage config file is required.")

        with open(sage_config_file) as f:
            param_data = json.load(f)

            # loading and set default parameters
            if "fragment_min_mz" in param_data["database"]:
                fragment_min_mz = param_data["database"]["fragment_min_mz"]
            else:
                param_data["database"]["fragment_min_mz"] = 150
                fragment_min_mz = 150
            if "fragment_max_mz" in param_data["database"]:
                fragment_max_mz = param_data["database"]["fragment_min_mz"]
            else:
                param_data["database"]["fragment_max_mz"] = 2000
                fragment_max_mz = 2000
            if "deisotope" in param_data:
                deisotope_bool = param_data["deisotope"]
            else:
                param_data["deisotope"] = True
                deisotope_bool = True
            if "min_matched_peaks" in param_data:
                min_matched_peaks = param_data["min_matched_peaks"]
            else:
                param_data["min_matched_peaks"] = 6
                min_matched_peaks = 6
            if "min_isotope_err" in param_data:
                min_isotope_err = param_data["min_isotope_err"]
            else:
                param_data["min_isotope_err"] = -1
                min_isotope_err = -1
            if "max_isotope_err" in param_data:
                max_isotope_err = param_data["max_isotope_err"]
            else:
                param_data["max_isotope_err"] = -1
                max_isotope_err = 3
            if "min_precursor_charge" in param_data:
                min_precursor_charge = param_data["min_precursor_charge"]
            else:
                min_precursor_charge = 2
                param_data["min_precursor_charge"] = 2
            if "max_precursor_charge" in param_data:
                max_precursor_charge = param_data["max_precursor_charge"]
            else:
                max_precursor_charge = 4
                param_data["max_precursor_charge"] = 4
            if "chimera" in param_data:
                chimera = param_data["chimera"]
            else:
                param_data["chimera"] = False
                chimera = False
            if "report_psms" in param_data:
                report_psms = param_data["report_psms"]
            else:
                param_data["report_psms"] = 1
                report_psms = 1
            if "annotate_matches" in param_data:
                annotate_matches = param_data["annotate_matches"]
            else:
                param_data["annotate_matches"] = False
                annotate_matches = False
            if "max_fragment_charge" in param_data:
                max_fragment_charge = param_data["max_fragment_charge"]
            else:
                param_data["max_fragment_charge"] = 1
                max_fragment_charge = 1
            if "search_thread" in param_data:
                search_thread = param_data["search_thread"]
            else:
                param_data["search_thread"] = 4
                search_thread = 4

        if fragment_tolerance != 0 and precursor_tolerance != 0 or not use_file_values:
            param_data["precursor_tol"]["ppm"] = [int(-1 * precursor_tolerance), int(precursor_tolerance)]
            precursor_tolerance = Tolerance(ppm=(int(-1 * precursor_tolerance), int(precursor_tolerance)))
            if fragment_type == "ppm":
                param_data["fragment_tol"][fragment_type] = [int(-1 * fragment_tolerance), int(fragment_tolerance)]
                fragment_tolerance = Tolerance(ppm=(int(-1 * fragment_tolerance), int(fragment_tolerance)))
            else:
                param_data["fragment_tol"][fragment_type] = [-1 * fragment_tolerance, fragment_tolerance]
                fragment_tolerance = Tolerance(da=(-1 * fragment_tolerance, fragment_tolerance))

        else:
            logging.info("Using the values from the file.")
            if "ppm" in param_data["precursor_tol"]:
                precursor_tolerance = param_data["precursor_tol"]["ppm"][1]
            else:
                precursor_tolerance = param_data["precursor_tol"]["da"][1]

            if "ppm" in param_data["fragment_tol"]:
                fragment_tolerance = param_data["fragment_tol"]["ppm"][1]
            else:
                fragment_tolerance = param_data["fragment_tol"]["da"][1]

        param_data["mzml_paths"] = mzml_files

        temp_sage_file = str(uuid.uuid4()) + ".json"
        with open(temp_sage_file, "w") as f:
            json.dump(param_data, f, indent=4)

        logging.info("Running SAGE with fragment tolerance: {} and precursor tolerance: {}".format(fragment_tolerance,
                                                                                                   precursor_tolerance))
        if not skip_spectra_preprocess:
            name_id_map = dict(zip(mzml_files, range(len(mzml_files))))
            with cf.ThreadPoolExecutor(search_thread) as executor:
                for m in mzml_files:
                    executor.submit(self.parse_mzml, m, name_id_map[m], param_data["database"]["fragment_min_mz"],
                                    param_data["database"]["fragment_max_mz"], deisotope_bool)

            # for future in cf.as_completed(tasks):
            #     all_queries.append(future.result())

        static_mods = param_data["database"]["static_mods"]
        variable_mods = param_data["database"]["variable_mods"]

        scorer = Scorer(report_psms=report_psms, min_matched_peaks=min_matched_peaks,
                        precursor_tolerance=precursor_tolerance,
                        fragment_tolerance=fragment_tolerance, min_isotope_err=min_isotope_err,
                        max_isotope_err=max_isotope_err,
                        min_precursor_charge=min_precursor_charge, max_precursor_charge=max_precursor_charge,
                        min_fragment_mass=fragment_min_mz, max_fragment_mass=fragment_max_mz, chimera=chimera,
                        annotate_matches=annotate_matches, max_fragment_charge=max_fragment_charge,
                        static_mods=static_mods, variable_mods=variable_mods
                        )

        indexed_db = self.generate_indexed_db(fasta_path=fasta_file, param_data=param_data)
        results = scorer.score_collection_psm(db=indexed_db, spectrum_collection=self.queries,
                                              num_threads=search_thread)

        psm_list = []
        for _, values in results.items():
            psm_list.extend(values)

        # after intensity a
        psm_list_rescore = rescore_lda(psm_list, verbose=False)
        psm_list_df = self.peptide_spectrum_match_list_to_pandas(psm_list_rescore, indexed_db, re_score=True)

        # TARGET = psm_list_df[psm_list_df.decoy == False]
        # DECOY = psm_list_df[psm_list_df.decoy]

        # plt.figure(figsize=(6, 4), dpi=150)
        #
        # plt.hist(TARGET.hyper_score, bins="auto", label="target")
        # plt.hist(DECOY.hyper_score, bins="auto", alpha=.5, label="decoy")
        # plt.legend()
        # plt.title("Target and decoy score distributions")
        # plt.xlabel("Score")
        # plt.ylabel("Count")
        # plt.show()

        TDC = target_decoy_competition_pandas(
            df=psm_list_df,
            # the method for TDC can be set here, e.g., PSM level, peptide level, or double competition (PSM and peptide level)
            method="peptide_psm_peptide",
            score="re_score"
        )

        # add PEP to the PSMs
        TDC["pep"] = calculate_pep(
            scores=TDC.score.values,
            decoys=TDC.decoy.values,
        )

        # after target decoy competition, hits can be filtered to, e.g., 1 percent expected FDR == q-value <= 0.01
        # TDC_filtered = TDC[(TDC.q_value <= 0.01) & (TDC.decoy == False)]
        RESULT = pd.merge(psm_list_df.drop(columns=["q_value", "score"]), TDC,
                          on=["spec_idx", "decoy", "match_idx"])

        # RESULT.insert(loc=0, column="file_name", value=os.path.splitext(os.path.basename(m))[0])
        RESULT.rename(columns={"q_value": "spectrum_q", "match_idx": "peptide"}, inplace=True)

        sage_table = pd.DataFrame(RESULT, index=None)

        sage_table = compute_entrapment_qvalues(sage_table)
        return sage_table


def extract_ptms(sage_table_target):
    def extract_modifications(peptide):
        return re.findall(r'([A-Z]\[\+[0-9.]+\])', peptide)

    # Apply the function to the peptides column and flatten the list of lists
    modifications = [mod for peptide in sage_table_target['peptide'] for mod in extract_modifications(peptide)]

    # Count the occurrences of each modification
    modification_counts = Counter(modifications)

    # convert modification counts to a dictionary key is the mod name, value is the count
    modification_counts_dict = dict(modification_counts)

    return modification_counts_dict


def get_stats_from_sage(sage_table: pd.DataFrame, fragment_tolerance: int = 0, precursor_tolerance: int = 0,
                        number_psms: int = 0) -> dict:
    sage_table = sage_table[sage_table['spectrum_q'] <= 0.01]
    sage_table = sage_table[sage_table['entrapment_qvalue'] <= 0.01]
    decoy_filter = sage_table['proteins'].str.contains("DECOY_")
    sage_table_target = sage_table[~decoy_filter]
    entrap_peptides = sage_table_target[sage_table_target['proteins'].str.contains("ENTRAP")]
    sage_table_decoy = sage_table[decoy_filter]
    precursor_std_error = sage_table_target['precursor_ppm'].std() * 4
    fragment_std_error = sage_table_target['fragment_ppm'].std() * 4

    plot_distribution(sage_table_target, sage_table_decoy, entrap_peptides, fragment_tolerance, precursor_tolerance)

    ptms_dist = extract_ptms(sage_table_target)

    stats = {
        "fragment_tolerance": fragment_tolerance,
        "precursor_tolerance": precursor_tolerance,
        "total_peptides": len(sage_table['peptide'].unique()),
        "total_entrap_peptides": len(entrap_peptides['peptide'].unique()),
        "total_decoys": len(sage_table_decoy['peptide'].unique()),
        "total_target": len(sage_table_target['peptide'].unique()),
        "precursor_std_error_plus4": precursor_std_error,
        "fragment_std_error_plus4": fragment_std_error,
        "number_psms": number_psms,
    }

    for ky, val in ptms_dist.items():
        stats[f"{ky}"] = val
        stats[f"{ky}_percentage"] = (val / len(sage_table_target['peptide'])) * 100

    return stats


def plot_distribution(sage_table_target, sage_table_decoy, entrap_peptides, fragment_tolerance, precursor_tolerance):
    plt.hist(sage_table_target['spectrum_q'], bins=20, alpha=0.5, label='Target')
    plt.hist(sage_table_decoy['spectrum_q'], bins=20, alpha=0.5, label='Decoy')
    plt.hist(entrap_peptides['spectrum_q'], bins=20, alpha=0.5, label='Entrap')
    number_psms = len(sage_table_target['peptide'])
    number_peptides = len(sage_table_target['peptide'].unique())
    plt.title("Q value distribution - Number of peptides: {} - Number of PSMs: {}".format(number_peptides, number_psms))
    plt.legend(loc='upper right')
    plt.xlabel(
        "Q value - fragment tolerance: {} - precursor tolerance: {}".format(fragment_tolerance, precursor_tolerance))
    file_name = folder_plots_uui + folder_sperator + "q_value_distribution-{}-{}.png".format(fragment_tolerance,
                                                                                             precursor_tolerance)
    plt.savefig(file_name)
    logging.info("Showing Q value distribution plot.")
    plt.close()
    # plt.show()


def compute_entrapment_qvalues(sage_data_frame: pd.DataFrame) -> pd.DataFrame:
    """ compute the q values using the entrapment peptides instead of decoys.
    :param sage_data_frame: panda data frame with the sage results
    """

    sage_data_frame.sort_values(by='re_score', ascending=False, inplace=True)

    # Use a more descriptive name for columns
    sage_data_frame['is_entrap'] = sage_data_frame['proteins'].apply(
        lambda x: 'entrap' if x.startswith('ENTRAP_') else 'target')
    sage_data_frame['is_decoy'] = sage_data_frame['proteins'].apply(
        lambda x: 'decoy' if x.startswith('DECOY_') else 'target')

    # Compute the q values using entrapment peptides
    sage_data_frame['spectrum_entrapment_q'] = 0
    sage_data_frame['cumulative_target'] = (sage_data_frame['is_entrap'] == 'target').cumsum()
    sage_data_frame['cumulative_entrap'] = (sage_data_frame['is_entrap'] == 'entrap').cumsum()
    sage_data_frame['FDR_ENTRAP'] = sage_data_frame['cumulative_entrap'] / sage_data_frame['cumulative_target']

    # Initialize the q-values with a large number
    sage_data_frame['entrapment_qvalue'] = 0

    # Use vectorized operation for calculating q-values
    sage_data_frame['entrapment_qvalue'] = sage_data_frame['FDR_ENTRAP'].expanding().min()

    return sage_data_frame


def compute_best_combination(sage_table: pd.DataFrame) -> int:
    """
    Compute the difference between the decoys and entrapmet qvalues.
    Sum all the differences and return the value.
    :param sage_table: pandas data frame with the sage results
    """
    # sage_data_frame['q_diff'] = (sage_data_frame['spectrum_q'] - sage_data_frame['spectrum_entrapment_q']).abs()
    # make copy of the data frame
    current_data_frame = sage_table.copy()
    current_data_frame = current_data_frame[current_data_frame['spectrum_q'] <= 0.01]
    current_data_frame = current_data_frame[current_data_frame['entrapment_qvalue'] <= 0.01]
    filter_decoy = current_data_frame['proteins'].str.contains("DECOY_")
    sage_table_target = current_data_frame[~filter_decoy]
    return len(sage_table_target['peptide'])


def combined_search(sagerunner, results: list = [], start_fragment_tolerance: int = 0,
                    start_precursor_tolerance: int = 0,
                    min_fragment_tolerance: int = 1, max_fragment_tolerance: int = 50,
                    min_precursor_tolerance: int = 10,
                    max_precursor_tolerance: int = 100, num_psms: int = 0, fragment_type: str = "ppm",
                    mzml_files: list = [],
                    grid_frag_steps: int = 10,
                    grid_prec_steps: int = 5, max_iterations=10, search_radius=1, fasta_file: str = None,
                    initial_temp=100,
                    cooling_rate=0.95, sage_config_file=None) -> Tuple[int, int, List[Dict]]:
    def acceptance_probability(old_value, new_value, temperature):
        if new_value > old_value:
            return 1.0
        return np.exp((new_value - old_value) / temperature)

    best_fragment_tolerance = int(start_fragment_tolerance)
    best_precursor_tolerance = int(start_precursor_tolerance)
    best_value = int(num_psms)

    # Coarse Grid Search
    precursor_tolerances = np.linspace(min_precursor_tolerance, max_precursor_tolerance, grid_prec_steps).astype(int)

    grid_best_value = 1
    for ft in precursor_tolerances:
        sage_table = sagerunner.run_sage(fragment_tolerance=start_fragment_tolerance, precursor_tolerance=ft,
                                         fragment_type=fragment_type, mzml_files=mzml_files, fasta_file=fasta_file,
                                         sage_config_file=sage_config_file, use_file_values=False,
                                         skip_spectra_preprocess=True)
        new_value = compute_best_combination(sage_table)
        results.append(get_stats_from_sage(sage_table, start_fragment_tolerance, ft, new_value))

        if (new_value - grid_best_value) / grid_best_value > 0.01:  # 1% improvement
            best_fragment_tolerance = start_fragment_tolerance
            best_precursor_tolerance = ft
            grid_best_value = new_value
            logging.info(
                "New Best value for precursor tolerance {}: {}".format(best_precursor_tolerance, grid_best_value))
        else:
            logging.info("Current value for precursor tolerance worst than previous {}: {}".format(ft, grid_best_value))
            break

    if best_precursor_tolerance == max_precursor_tolerance:
        logging.info(
            "Best precursor tolerance is the maximum value {}. You have to consider locking for for unknown modifications.".format(
                best_precursor_tolerance))

    fragment_tolerances = np.linspace(start_fragment_tolerance, max_fragment_tolerance, grid_frag_steps).astype(int)

    for ft in fragment_tolerances:
        for pt in precursor_tolerances:
            sage_table = sagerunner.run_sage(fragment_tolerance=ft, precursor_tolerance=pt,
                                             fragment_type=fragment_type, mzml_files=mzml_files, fasta_file=fasta_file,
                                             sage_config_file=sage_config_file, skip_spectra_preprocess=True,
                                             use_file_values=False)
            new_value = compute_best_combination(sage_table)
            results.append(get_stats_from_sage(sage_table, ft, pt, new_value))

            if (new_value - grid_best_value) / grid_best_value > 0.01:  # 1% improvement
                best_fragment_tolerance = ft
                best_precursor_tolerance = pt
                best_value = new_value
                logging.info("New Best value for fragment tolerance {}: {}".format(best_fragment_tolerance, best_value))
            else:
                logging.info(
                    "Current value for fragment tolerance worst than previous {}: {}".format(ft, grid_best_value))
                break

    # Simulated Annealing
    if grid_best_value < best_value:
        best_value = int(num_psms)
        logging.info("Best value from grid search is better than the initial value. Using it.")
        best_fragment_tolerance = start_fragment_tolerance
        best_precursor_tolerance = start_precursor_tolerance

    temperature = initial_temp
    current_fragment_tolerance = best_fragment_tolerance
    current_precursor_tolerance = best_precursor_tolerance

    for _ in range(max_iterations):
        fragment_tolerance = np.random.randint(max(min_fragment_tolerance, current_fragment_tolerance - search_radius),
                                               min(max_fragment_tolerance, current_fragment_tolerance + search_radius))
        precursor_tolerance = np.random.randint(
            max(min_precursor_tolerance, current_precursor_tolerance - search_radius),
            min(max_precursor_tolerance, current_precursor_tolerance + search_radius))

        sage_table = sagerunner.run_sage(fragment_tolerance=fragment_tolerance, precursor_tolerance=precursor_tolerance,
                                         fragment_type=fragment_type, mzml_files=mzml_files, fasta_file=fasta_file,
                                         sage_config_file=sage_config_file, skip_spectra_preprocess=True,
                                         use_file_values=False)
        new_value = compute_best_combination(sage_table)
        results.append(get_stats_from_sage(sage_table, fragment_tolerance, precursor_tolerance, new_value))

        if acceptance_probability(best_value, new_value, temperature) > random.random():
            current_fragment_tolerance = fragment_tolerance
            current_precursor_tolerance = precursor_tolerance
            best_fragment_tolerance = fragment_tolerance
            best_precursor_tolerance = precursor_tolerance
            best_value = new_value

        temperature *= cooling_rate

    return best_fragment_tolerance, best_precursor_tolerance, results


@click.command("tolerances", help="Optimize the fragment and precursor tolerances using the SAGE algorithm.")
@click.option("--fragment-type", default="ppm", help="The type of fragment tolerance to use (ppm or da)")
@click.option("--mzml-path", default=".", help="The path to the mzML files to use for the SAGE analysis.")
@click.option("--initial-fragment-tolerance", help="The initial fragment tolerance to use for the optimization.",
              type=int, default=20)
@click.option("--initial-precursor-tolerance", help="The initial precursor tolerance to use for the optimization.",
              type=int, default=20)
@click.option("--min-fragment-tolerance", default=1, help="The minimum fragment tolerance to consider.", type=int)
@click.option("--max-fragment-tolerance", default=50, help="The minimum precursor tolerance to consider.", type=int)
@click.option("--min-precursor-tolerance", default=10, help="The minimum precursor tolerance to consider.", type=int)
@click.option("--max-precursor-tolerance", default=50, help="The maximum precursor tolerance to consider.", type=int)
@click.option("--fasta-file", default="Homo-sapiens-uniprot-reviewed-contaminants-entrap-decoy-20240615.fasta",
              help="The path to the fasta file to use for the SAGE analysis.")
@click.option("--sage-config-file", default="general-sage.json",
              help="The path to the Sage config file to use for the SAGE analysis.")
@click.option("--max-iterations", default=10, help="The maximum number of iterations to run the optimization.",
              type=int)
def tolerances(fragment_type: str, mzml_path: str, initial_fragment_tolerance: int, initial_precursor_tolerance: int,
               min_fragment_tolerance: int, max_fragment_tolerance: int, min_precursor_tolerance: int,
               max_precursor_tolerance: int,
               fasta_file, sage_config_file, max_iterations: int):
    results = []

    # detect absolute all the mzML files in the mzml-path
    mzml_files = [os.path.join(mzml_path, f) for f in os.listdir(mzml_path) if f.endswith(".mzML")]

    # generate uui unique identifier for the folder with the plots
    global folder_plots_uui
    folder_plots_uui = str(uuid.uuid4())
    os.makedirs(folder_plots_uui)

    num_psms = 0
    best_precursor_tolerance = 0
    best_fragment_tolerance = 0

    if sage_config_file is None:
        sage_config_file = "general-sage.json"

    if initial_fragment_tolerance is not None and initial_precursor_tolerance is not None:

        if os.path.exists(fasta_file):
            sagerunner = SageRunner()
        else:
            logging.error(f"File {fasta_file} does not exist.")
            raise FileNotFoundError(f"File {fasta_file} does not exist.")

        sage_table = sagerunner.run_sage(int(initial_fragment_tolerance), int(initial_precursor_tolerance),
                                         fragment_type,
                                         mzml_files, sage_config_file=sage_config_file, fasta_file=fasta_file,
                                         use_file_values=False)
        sage_table.to_csv("sage_table.csv", index=False)
        num_psms = compute_best_combination(sage_table)
        stats = get_stats_from_sage(sage_table, initial_fragment_tolerance, initial_precursor_tolerance, num_psms)
        results.append(stats)
        initial_fragment_tolerance = math.ceil(stats["fragment_std_error_plus4"])
        best_precursor_tolerance = initial_precursor_tolerance
        best_fragment_tolerance = initial_fragment_tolerance

    best_fragment_tolerance, best_precursor_tolerance, results = combined_search(sagerunner, results=results,
                                                                                 start_fragment_tolerance=best_fragment_tolerance,
                                                                                 start_precursor_tolerance=best_precursor_tolerance,
                                                                                 min_fragment_tolerance=initial_fragment_tolerance,
                                                                                 max_fragment_tolerance=max_fragment_tolerance,
                                                                                 min_precursor_tolerance=min_precursor_tolerance,
                                                                                 max_precursor_tolerance=max_precursor_tolerance,
                                                                                 num_psms=num_psms,
                                                                                 fragment_type=fragment_type,
                                                                                 mzml_files=mzml_files,
                                                                                 max_iterations=max_iterations,
                                                                                 fasta_file=fasta_file,
                                                                                 sage_config_file=sage_config_file)

    print("Best tolerances found: Fragment tolerance: {} - Precursor tolerance: {}".format(best_fragment_tolerance,
                                                                                           best_precursor_tolerance))

    # Write the results to a file
    results_df = pd.DataFrame(results)

    # Aggregate the data to handle duplicate entries by averaging number_psms
    df_agg = results_df.groupby(['fragment_tolerance', 'precursor_tolerance']).agg(
        {'number_psms': 'mean'}).reset_index()

    # Pivot the aggregated dataframe to prepare for plotting
    df_pivot = df_agg.pivot(index='fragment_tolerance', columns='precursor_tolerance', values='number_psms')

    # Plot the results
    ax = df_pivot.plot(kind='line', marker='o', figsize=(10, 6))

    # Set plot title and labels
    ax.set_title('Number of PSMs vs. Fragment Tolerances')
    ax.set_xlabel('Fragment Tolerance')
    ax.set_ylabel('Number of PSMs')

    # Set legend title
    ax.legend(title='Precursor Tolerances')

    # Save the plot to a file
    output_file = f"{folder_plots_uui}/final_results_tolerances.png"
    plt.savefig(output_file)
    plt.close()

    print(f"Plot saved to {output_file}")

    results_df.to_csv("sage_results_tolerances.tsv", sep="\t", index=False)


@click.command("ptms", help="Extract PTMs from SAGE results.")
@click.option("--mzml-path", default=".", help="The path to the mzML files to use for the SAGE analysis.")
@click.option("--fasta-file", default="Homo-sapiens-uniprot-reviewed-contaminants-entrap-decoy-20240615.fasta",
              help="The path to the fasta file to use for the SAGE analysis.")
@click.option("--sage-config-file", default="general-sage-ptms.json",
              help="The path to the Sage config file to use for the SAGE analysis.", required=True)
def ptms(mzml_path: str, fasta_file: str, sage_config_file: str):
    # detect absolute all the mzML files in the mzml-path

    if not os.path.exists(mzml_path):
        logging.error(f"Folder {mzml_path} does not exist.")
        raise FileNotFoundError(f"Folder {mzml_path} does not exist.")

    if not os.path.exists(fasta_file):
        logging.error(f"File {fasta_file} does not exist.")
        raise FileNotFoundError(f"File {fasta_file} does not exist.")

    if not os.path.exists(sage_config_file):
        logging.error(f"File {sage_config_file} does not exist.")
        raise FileNotFoundError(f"File {sage_config_file} does not exist.")
    else:
        sage_params = json.load(open(sage_config_file))
        if "ppm" in sage_params["precursor_tol"]:
            precursor_tolerance = sage_params["precursor_tol"]["ppm"][1]
        else:
            precursor_tolerance = sage_params["precursor_tol"]["da"][1]

        if "ppm" in sage_params["fragment_tol"]:
            fragment_tolerance = sage_params["fragment_tol"]["ppm"][1]
        else:
            fragment_tolerance = sage_params["fragment_tol"]["da"][1]

    mzml_files = [os.path.join(mzml_path, f) for f in os.listdir(mzml_path) if f.endswith(".mzML")]

    global folder_plots_uui
    folder_plots_uui = str(uuid.uuid4())
    os.makedirs(folder_plots_uui)
    sagerunner = SageRunner()
    sage_table = sagerunner.run_sage(fragment_type="ppm", mzml_files=mzml_files, fasta_file=fasta_file,
                                     sage_config_file=sage_config_file, use_file_values=True)

    num_psms = compute_best_combination(sage_table)
    stats = get_stats_from_sage(sage_table, fragment_tolerance=fragment_tolerance,
                                precursor_tolerance=precursor_tolerance, number_psms=num_psms)

    # check verbose the PTMs that percentace bigger than 1% of the psms
    logging.info("PTMs with percentage bigger than 1% of the psms:")
    for key, value in stats.items():
        if "_percentage" in key and value > 1:
            logging.info(f"{key} - {value}")

    # print dictionary stats to a file as key value.
    with open("sage_stats_ptms.tsv", "w") as f:
        for key, value in stats.items():
            f.write(f"{key}\t{value}\n")


def sage2PTMShepherd(sage_table):
    ptmsherd_schema = ["Spectrum", "Spectrum File", "Peptide", "Modified Peptide", "Peptide Length", "Charge",
                       "Retention", "Observed Mass", "Calibrated Observed Mass",
                       "Observed M/Z", "Calibrated Observed M/Z", "Calculated Peptide Mass", "Calculated M/Z",
                       "Delta Mass", "Expectation", "Hyperscore", "Nextscore", "Probability",
                       "Number of Missed Cleavages", "Intensity",
                       "Assigned Modifications", "Observed Modifications",
                       "Purity", "Protein", "Protein ID"]
    psm_df = pd.DataFrame(columns=ptmsherd_schema)
    unimod_pattern = re.compile(r"\[UNIMOD:\d+]")
    db = mass.Unimod()
    new_rows = []
    for _, row in sage_table.iterrows():
        Spectrum = row["spec_idx"]
        SpectrumFile = row["spec_idx"].split(".")[0]
        peptide = unimod_pattern.sub(repl="", string=row['peptide'])
        modifications = unimod_pattern.finditer(row["peptide"])
        modified_peptide = row["peptide"]
        AssignedModifications = []
        previous_idx = 0
        for modification in modifications:
            position = modification.span()[0] - 1 - previous_idx
            unimod_term = modification.group(0).replace("[UNIMOD:", "").replace("]", "")
            mass_da = db.by_id(unimod_term)["mono_mass"]
            AssignedModifications.append(str(position + 1) + str(peptide[position]) + "(" + str(mass_da) + ")")
            MassResiduePlusMod = int(mass.fast_mass(peptide[position]) + mass_da - 18.001)
            modified_peptide = modified_peptide.replace(peptide[position] + modification.group(0),
                                                        peptide[position] + "[" + str(MassResiduePlusMod) + "]")
            previous_idx += len(modification.group(0))
        if modified_peptide == peptide:
            modified_peptide = ""
        AssignedModifications = ", ".join(AssignedModifications)

        PeptideLength = len(peptide)
        Charge = row["charge"]
        Retention = row["retention_time_observed"] * 60
        ObservedMass = row["mono_mass_observed"]
        ObservedMZ = row["mono_mass_observed"] / row["charge"] + 1.007276
        CalibratedObservedMZ = ObservedMZ
        CalculatedPeptideMass = row["mono_mass_calculated"]
        CalibratedObservedMass = ObservedMass
        CalculatedMZ = row["mono_mz_calculated"]
        DeltaMass = row["delta_mass"]
        Hyperscore = row["hyperscore"]
        Probability = 1 - row["spectrum_q"]
        Nextscore = row["delta_next"]
        NumberofMissedCleavages = row["missed_cleavages"]
        Intensity = row["intensity_ms1"]
        Protein = row["proteins"]
        ProteinID = row["proteins"].split("|")[1]
        new_rows.append({"Spectrum": Spectrum, "Spectrum File": SpectrumFile,
                         "Peptide": peptide,
                         "Modified Peptide": modified_peptide,
                         "Peptide Length": PeptideLength, "Charge": Charge,
                         "Retention": Retention, "Observed Mass": ObservedMass,
                         "Calibrated Observed Mass": CalibratedObservedMass,
                         "Observed M/Z": ObservedMZ,
                         "Calibrated Observed M/Z": CalibratedObservedMZ,
                         "Calculated Peptide Mass": CalculatedPeptideMass,
                         "Calculated M/Z": CalculatedMZ, "Delta Mass": DeltaMass,
                         "Hyperscore": Hyperscore, "Nextscore": Nextscore,
                         "Probability": Probability,
                         "Number of Missed Cleavages": NumberofMissedCleavages,
                         "Intensity": Intensity,
                         "Assigned Modifications": AssignedModifications,
                         "Protein": Protein, "Protein ID": ProteinID})

    psm_df = pd.concat([psm_df, pd.DataFrame.from_records(new_rows)],
                       ignore_index=True)

    psm_df.to_csv("inputPTMShepherd.tsv", index=False, sep="\t")


@click.group(context_settings=CONTEXT_SETTINGS)
def cli():
    """
    Main function to run the CLI
    """
    pass


cli.add_command(tolerances)
cli.add_command(ptms)

if __name__ == "__main__":
    cli()
