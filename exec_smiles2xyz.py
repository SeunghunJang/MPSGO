from utils import process_excel_to_xyz, create_output_folder

# Define paths for the Excel file, output folder, and error log file
excel_file_path = "Data/Sample.xlsx"
output_folder_path = "Gen_XYZ_files"
error_file_path = "Log_error.txt"

# Create the folder if it doesn't exist
create_output_folder(output_folder_path)

# Process the Excel file to generate XYZ files and log any errors
process_excel_to_xyz(excel_file_path, output_folder_path, error_file_path)
