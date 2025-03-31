import os
import re
import csv
from multiprocessing import Pool

# Set directory path.
base_directory = 'Gen_molecular_pairs_cal'
output_filepath = 'Gen_molecular_pairs_analysis/Properties.csv'

# Extracts energy values and related information.
def extract_energy_values(filepath):
    total_energy = None
    homo_lumo_gap = None
    total_run_time = None
    homo = None
    lumo = None
    gradient_norm = None

    with open(filepath, 'r') as file:
        for line in file:
            # Search for total energy in Hartree
            if 'TOTAL ENERGY' in line:
                match = re.search(r"TOTAL ENERGY\s+([-+]?\d*\.\d+)\s*Eh", line)
                if match:
                    total_energy = float(match.group(1))

            # Search for HOMO-LUMO gap in eV
            elif 'HOMO-LUMO GAP' in line:
                match = re.search(r"HOMO-LUMO GAP\s+([-+]?\d*\.\d+)\s*eV", line)
                if match:
                    homo_lumo_gap = float(match.group(1))

            # Search for total run time and convert it to seconds
            elif 'TOTAL RUN TIME' in line:
                match = re.search(r"TOTAL RUN TIME:\s+(\d+ days \d+ hours \d+ minutes \d+ seconds \d+ msec)", line)
                if match:
                    total_run_time = convert_time_to_seconds(match.group(1))

            # Search for gradient norm in Hartree per Bohr radius
            elif 'GRADIENT NORM' in line:
                match = re.search(r"GRADIENT NORM\s+([-+]?\d*\.\d+)\s*Eh/α", line)
                if match:
                    gradient_norm = float(match.group(1))

            # Search for HOMO and LUMO values
            else:
                homo_match = re.search(r"^\s*\d+\s+\S+\s+([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)\s+\(HOMO\)", line)
                lumo_match = re.search(r"^\s*\d+\s+\S*\s+([-+]?\d*\.\d+)\s+([-+]?\d*\.\d+)\s+\(LUMO\)", line)
                if homo_match:
                    homo = float(homo_match.group(2))
                if lumo_match:
                    lumo = float(lumo_match.group(2))

    return total_energy, homo_lumo_gap, total_run_time, gradient_norm, homo, lumo

# Converts a time string into total seconds
def convert_time_to_seconds(time_str):
    match = re.search(r"(\d+) days (\d+) hours (\d+) minutes (\d+) seconds (\d+) msec", time_str)
    if match:
        days = int(match.group(1))
        hours = int(match.group(2))
        minutes = int(match.group(3))
        seconds = int(match.group(4))
        msec = int(match.group(5))
        total_seconds = days * 86400 + hours * 3600 + minutes * 60 + seconds + msec / 1000.0
        return total_seconds
    return None

# Extracts energy and runtime data from all std.out files in a folder
def extract_from_folder(folder):
    print(f"Processing folder: {folder}")
    results = []
    for root, dirs, files in os.walk(folder):
        for file in files:
            if file == 'std.out':
                filepath = os.path.join(root, file)
                total_energy, homo_lumo_gap, total_run_time, gradient_norm, homo, lumo = extract_energy_values(filepath)

                # Only add the result if any value was found
                if total_energy is not None or homo_lumo_gap is not None or total_run_time is not None or gradient_norm is not None or homo is not None or lumo is not None:
                    results.append({
                        'folder': root,
                        'total_energy': total_energy,
                        'homo_lumo_gap': homo_lumo_gap,
                        'homo': homo,
                        'lumo': lumo,
                        'gradient_norm': gradient_norm,
                        'total_run_time': total_run_time
                    })

    print(f"Finished processing folder: {folder}")
    return results


# Extract data from all folders with multiprocessing
def extract_from_all_folders(base_directory):
    # Get list of folders to process
    folders = [os.path.join(base_directory, folder) for folder in os.listdir(base_directory) if os.path.isdir(os.path.join(base_directory, folder))]

    print(f"Total folders to process: {len(folders)}")

    # Use multiprocessing to process folders in parallel
    with Pool() as pool:
        results = pool.map(extract_from_folder, folders)

    # Flatten list of results from all folders
    return [item for sublist in results for item in sublist]


# Writes the extracted results to a CSV.
def write_results_to_csv(results, output_filepath):
    os.makedirs(os.path.dirname(output_filepath), exist_ok=True)

    print(f"Writing results to {output_filepath}...")
    with open(output_filepath, 'w', newline='') as csvfile:
        fieldnames = ['Folder Name', 'TOTAL ENERGY [Eh]', 'HOMO-LUMO GAP [eV]', 'HOMO [eV]', 'LUMO [eV]', 'GRADIENT NORM [Eh/a0]', 'TOTAL RUN TIME [sec]']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for result in results:
            writer.writerow({
                'Folder Name': result['folder'],
                'TOTAL ENERGY [Eh]': result['total_energy'],
                'HOMO-LUMO GAP [eV]': result['homo_lumo_gap'],
                'HOMO [eV]': result['homo'],
                'LUMO [eV]': result['lumo'],
                'GRADIENT NORM [Eh/a0]': result['gradient_norm'],
                'TOTAL RUN TIME [sec]': result['total_run_time']
            })

    print("Results have been written to CSV.")


# Main execution
if __name__ == "__main__":
    # Extract data
    results = extract_from_all_folders(base_directory)

    # Write results to CSV
    write_results_to_csv(results, output_filepath)

    print(f"Process completed. Results have been written to {output_filepath}")
