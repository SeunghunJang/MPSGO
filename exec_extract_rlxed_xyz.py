import os
import shutil

# Define directories
source_dir = './Gen_molecular_pairs_cal'
destination_dir = './Gen_molecular_pairs_xyz'

# Create directory if it doesn't exist
os.makedirs(destination_dir, exist_ok=True)

# List subdirectories in source dir
subdirs = os.listdir(source_dir)
num_subdirs = len(subdirs)
print(f"Total subdirectories to process: {num_subdirs}")

# Loop through each subdirectory
for subdir in subdirs:
    print(f"Processing: {subdir}")
    
    # Paths for input.xyz file
    input_xyz_path = os.path.join(source_dir, subdir, 'input.xyz')

    # Use the subdirectory name as the new file name
    new_file_name = f"{subdir.replace('mol_', '')}"
    subdir_dest_path = os.path.join(destination_dir, new_file_name)
    
    # Check if input.xyz exists and copy it
    if os.path.isfile(input_xyz_path):
        try:
            shutil.copy(input_xyz_path, subdir_dest_path)
            print(f"Copied {input_xyz_path} to {subdir_dest_path}")
        except Exception as e:
            print(f"Error copying {input_xyz_path}: {e}")
    else:
        print(f"File {input_xyz_path} not found.")

