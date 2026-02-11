import sys
import os
import subprocess
import click

@click.command("param-medic", short_help="Process RAW files with param-medic tool to predict precursor and fragment tolerances")
@click.option("--input-file", type=click.Path(exists=True), required = True)
def main(input_file: str):

    # Check if the input file exists
    if not os.path.isfile(input_file):
        print(f"Error: {input_file} is not a valid file")
        sys.exit(1)

    file_mzml = input_file
    if input_file.lower().endswith('.raw'):
        # Convert raw file to mzML
        print("Converting the raw file into mzML")
        subprocess.run(['ThermoRawFileParser.sh', f'-i={input_file}', '-o=./', '-f=2'], check=True)
        file_mzml = input_file.replace('.raw', '.mzML')    # Replace the psi term of the software msconvert with thermorawfileparser

    print("Replacing psi term for msconvert to ThermoRawFileParser")
    input_mzml = f"{os.path.splitext(file_mzml)[0]}.mzML"
    fixed_mzml = input_file.replace('.mzML', '_fixed.mzML')
    with open(input_mzml, 'r') as file:
        data = file.read()
    data = data.replace('1003145', '1000615')
    with open(fixed_mzml, 'w') as file:
        file.write(data)

    # Run param medic to guess precursor and fragment tolerances
    print("Running param medic to predict the precursor and fragment tolerances")
    subprocess.run(
        ['crux', 'param-medic', fixed_mzml, '--overwrite', 'T', '--pm-charges', '2,3,4', '--pm-min-peak-pairs', '20'],
        check=True)

    result_file_param_medic = "crux-output/param-medic.log.txt"
    precursor_error = None
    fragment_error = None
    if os.path.isfile(result_file_param_medic):
        with open(result_file_param_medic, 'r') as file:
            param = file.read()
            # Extract the precursor and fragment tolerances
            # example: INFO: Precursor error: 39.10 ppm
            # INFO: Fragment bin size: 0.0056 Th
            precursor_error = param.split('Precursor error: ')[1]
            fragment_error = param.split('Fragment bin size: ')[1]

    # Delete the intermediate files
    print("Deleting the intermediate files")
    os.remove(fixed_mzml)
    if input_file.lower().endswith('.raw'):
        os.remove(input_mzml)

    print(f"Modified file created: {os.path.splitext(input_file)[0]}_modified.{os.path.splitext(input_file)[1]}")

    if precursor_error and fragment_error:
        print(f"Precursor error: {precursor_error}")
        print(f"Fragment error: {fragment_error}")


if __name__ == "__main__":
    main()
