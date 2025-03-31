import os
from tqdm import tqdm
from utils import read_xyz_file, create_output_folder, merge_molecules

# Input and Output Folders
input_folder = "Gen_XYZ_files"
output_folder = "Gen_molecular_pairs"

# Create the output folder if it doesn't exist
create_output_folder(output_folder)

# Get list of files in the input folder
input_files = os.listdir(input_folder)

# Group files into molecule pairs based on naming conventions
molecule_pairs = {}
for file in input_files:
    if file.endswith("_1.xyz"):
        molecule_pairs.setdefault(file[:-6], []).append(file)
    elif file.endswith("_2.xyz"):
        molecule_pairs.setdefault(file[:-6], []).append(file)

# Merge molecules for each pair
for pair_name, pair_files in tqdm(molecule_pairs.items()):
    if len(pair_files) != 2:
        print(f"Error: Incomplete pair files for pair {pair_name}")
        continue

    file_path_1 = os.path.join(input_folder, pair_files[0])
    file_path_2 = os.path.join(input_folder, pair_files[1])

    symbols_1, molecule_coordinates_1 = read_xyz_file(file_path_1)
    symbols_2, molecule_coordinates_2 = read_xyz_file(file_path_2)

    merge_count = merge_molecules(pair_name, symbols_1, molecule_coordinates_1, symbols_2, molecule_coordinates_2, output_folder=output_folder)

    if merge_count < 10:
        print(f"Less than 10 merged XYZ files were created for pair {pair_name}.")
    else:
        print(f"Merged XYZ files created successfully for pair {pair_name}.")

print("Merge process completed.")
