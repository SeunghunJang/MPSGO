import os
import shutil
import subprocess

# Define directories and file path
input_dir = 'Gen_molecular_pairs'
input_path = './Input_dft/input.inp'
qsub_path = './Input_dft/qsub.sh'

# Count number of xyz files in input directory
files = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f))]
num_files = len(files)
print(f"Number of xyz files: {num_files}")

# List files in the input directory
mpt = os.listdir(input_dir)

# Loop through the specified range of files
for i in range(0, num_files):
    xyz_name = mpt[i]
    print(f"Processing {xyz_name}")

    # Create the directory for the current job
    job_dir = f'./{input_dir}_cal/{xyz_name}'
    os.makedirs(job_dir, exist_ok=True)

    # Define file paths
    xyz_path = f'./{input_dir}/{xyz_name}'

    # Copy the necessary files if they exist
    for src, dest in [(xyz_path, f'{job_dir}/mol.xyz'),
                      (input_path, f'{job_dir}/input.inp'),
                      (qsub_path, f'{job_dir}/qsub.sh')]:

        try:
            if os.path.exists(src):
                shutil.copy(src, dest)
            else:
                print(f"Error: {src} not found.")
                break
        except Exception as e:
            print(f"Error copying {src} to {dest}: {e}")
            break
    else:
        # Submit the job if all files were copied successfully
        os.chdir(job_dir)

        try:
            # Run the sbatch command and capture stdout and stderr
            result = subprocess.run(['sbatch', 'qsub.sh'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)

            # Check the result of the job submission
            if result.returncode != 0:
                print(f"Error submitting job: {result.stderr}")
            else:
                print(f"Job submitted successfully: {result.stdout}")

        except Exception as e:
            print(f"Error during job submission: {e}")
        finally:
            # Ensure we return to the parent directory
            os.chdir('../../')
