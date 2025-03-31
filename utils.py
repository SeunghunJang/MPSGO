import os
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

# Add hydrogen atoms to the molecule
def add_hydrogens(mol):
    mol_with_hydrogens = Chem.AddHs(mol)
    return mol_with_hydrogens

# Write the 3D coordinates of a molecule to an XYZ file.
def write_mol_xyz_file(mol, smiles, xyz_file_path):
    with open(xyz_file_path, 'w') as xyz_file:
        num_atoms = mol.GetNumAtoms()           # Get number of atoms in the molecule
        xyz_file.write(f"{num_atoms}\n")        # Write atom count in the first line
        xyz_file.write(f"{smiles}\n")  # Optional header line

        # Write atom details (symbol and coordinates)
        for atom in mol.GetAtoms():
            pos = mol.GetConformer().GetAtomPosition(atom.GetIdx())
            symbol = atom.GetSymbol()
            xyz_file.write(f"{symbol} {pos.x} {pos.y} {pos.z}\n")

# Convert two SMILES strings to XYZ files.
def smiles_to_xyz_file(smiles_1, smiles_2, output_folder, xyz_filename, error_file):
    try:
        # Convert SMILES to RDKit Mol objects
        mol_1 = Chem.MolFromSmiles(smiles_1)
        mol_2 = Chem.MolFromSmiles(smiles_2)

        # Check if both molecules were successfully parsed
        if mol_1 is not None and mol_2 is not None:
            # Add hydrogen atoms to both molecules
            mol_with_hydrogens_1 = add_hydrogens(mol_1)
            mol_with_hydrogens_2 = add_hydrogens(mol_2)

            # Generate 3D structures and optimize them using MMFF94 force field
            AllChem.EmbedMolecule(mol_with_hydrogens_1)
            AllChem.MMFFOptimizeMolecule(mol_with_hydrogens_1)
            AllChem.EmbedMolecule(mol_with_hydrogens_2)
            AllChem.MMFFOptimizeMolecule(mol_with_hydrogens_2)

            # Define output paths for XYZ files
            xyz_file_path_1 = f"{output_folder}/{xyz_filename}_1.xyz"
            xyz_file_path_2 = f"{output_folder}/{xyz_filename}_2.xyz"

            # Write the 3D coordinates to XYZ files
            write_mol_xyz_file(mol_with_hydrogens_1, smiles_1, xyz_file_path_1)
            write_mol_xyz_file(mol_with_hydrogens_2, smiles_2, xyz_file_path_2)

            print(f"XYZ files saved: {xyz_file_path_1}, {xyz_file_path_2}")
        else:
            print(f"Invalid SMILES: {smiles_1} or {smiles_2}")
    except Exception as e:
        print(f"Error processing SMILES: {smiles_1} or {smiles_2}")
        print(f"Error message: {str(e)}")

        # Log errors to the error file
        with open(error_file, 'a') as err_file:
            err_file.write(f"{smiles_1}, {smiles_2}\n")

# Reads an Excel file with SMILES pairs and generates XYZ files
def process_excel_to_xyz(excel_file, output_folder, error_file):
    df = pd.read_excel(excel_file)

    if 'SMILES_1' in df.columns and 'SMILES_2' in df.columns:
        # Iterate over each row in the DataFrame
        for index, row in df.iterrows():
            smiles_value_1 = row['SMILES_1']
            smiles_value_2 = row['SMILES_2']
            xyz_filename = f"mol_{index + 1}"  # Create a unique filename for each pair
            # Generate XYZ files for the SMILES pair
            smiles_to_xyz_file(smiles_value_1, smiles_value_2, output_folder, xyz_filename, error_file)
    else:
        print("SMILES_1 or SMILES_2 column not found in the Excel file.")

# Reads atomic coordinates from an XYZ file.
def read_xyz_file(file_path):
    symbols = []
    coordinates = []
    with open(file_path, 'r') as f:
        num_atoms = int(f.readline())        # Read the number of atoms
        f.readline()                         # Skip the comment line
        for _ in range(num_atoms):
            line = f.readline().split()
            symbol = line[0]                 # Atomic symbol
            x, y, z = map(float, line[1:4])  # Atomic coordinates
            symbols.append(symbol)
            coordinates.append([x, y, z])
    return symbols, np.array(coordinates)

# Writes atomic coordinates and a distance value to an XYZ file.
def write_xyz_file(file_path, symbols, coordinates, min_distance):
    with open(file_path, 'w') as f:
        f.write(f"{len(symbols)}\n")                                      # Number of atoms
        f.write(f"Minimum distance between molecules: {min_distance}\n")  # Distance information
        for symbol, coord in zip(symbols, coordinates):
            f.write(f"{symbol} {coord[0]} {coord[1]} {coord[2]}\n")       # Atomic symbol and coordinates

# Creates the output folder if it doesn't exist.
def create_output_folder(output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

# Finds the center and radius of the sphere enclosing the molecule.
def find_enclosing_sphere(coordinates):
    center = np.mean(coordinates, axis=0)                                # Center of the molecule
    distances = np.linalg.norm(coordinates - center, axis=1)             # Distances from center
    radius = np.max(distances)                                           # Radius is the maximum distance
    return center, radius

# Calculates the minimum distance between two sets of atomic coordinates.
def calculate_distance(coord1, coord2):
    min_distance = float('inf')                                          # Initialize with a very large value
    for atom1 in coord1:
        for atom2 in coord2:
            distance = np.linalg.norm(atom1 - atom2)                     # Euclidean distance between two atoms
            if distance < min_distance:
                min_distance = distance                                  # Update if a smaller distance is found
    return min_distance

# Generates random points on the surface of an ellipsoid based on the first molecule's coordinates.
def generate_random_points(coordinates_1, radius_2, num_points=20):
    def fit_ellipsoid_to_points(points):
        centroid = np.mean(points, axis=0)                             # Centroid of the points
        centered_points = points - centroid                            # Center the points
        covariance_matrix = np.cov(centered_points, rowvar=False)      # Covariance matrix
        eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)  # Eigen decomposition
        sort_indices = np.argsort(eigenvalues)[::-1]                   # Sort eigenvalues in descending order
        eigenvalues = eigenvalues[sort_indices]
        eigenvectors = eigenvectors[:, sort_indices]
        radii = np.sqrt(eigenvalues)                                   # Calculate radii as the square root of eigenvalues
        return centroid, eigenvectors, radii

    def generate_random_points_on_ellipsoid_surface(centroid, eigenvectors, radii, num_points=num_points, max_iter=1000,
                                                    radius_2=radius_2):
        points = []
        num_generated = 0
        while num_generated < num_points:
            phi = np.random.uniform(0, 2 * np.pi)                      # Azimuthal angle
            theta = np.random.uniform(0, np.pi)                        # Polar angle
            x = (radii[0] + radius_2) * np.cos(phi) * np.sin(theta) + centroid[0]
            y = (radii[1] + radius_2) * np.sin(phi) * np.sin(theta) + centroid[1]
            z = (radii[2] + radius_2) * np.cos(theta) + centroid[2]
            point = np.array([x, y, z])
            diff = point - centroid
            normalized_point = np.dot(diff, eigenvectors)              # Normalize the point using the eigenvectors
            normalized_point /= radii + radius_2
            distance = np.linalg.norm(normalized_point)                # Check if the point lies on the surface
            if np.abs(distance - 1) < 1e-2:
                points.append(point)
                num_generated += 1
            max_iter -= 1
            if max_iter <= 0:
                break
        return np.array(points)

    centroid, eigenvectors, radii = fit_ellipsoid_to_points(coordinates_1) # Fit ellipsoid to molecule 1
    max_atom_distance = np.max(np.linalg.norm(coordinates_1 - centroid, axis=1))  # Maximum distance of atoms from centroid
    radii += max_atom_distance                                             # Adjust radii based on the atom distances
    random_points = generate_random_points_on_ellipsoid_surface(centroid, eigenvectors, radii, radius_2=radius_2)
    return random_points

# Merges two molecules by placing one molecule at various positions relative to the other.
def merge_molecules(file_name, symbols_1, coordinates_1, symbols_2, coordinates_2, radius_2_offset=0.0, num_points=20, output_folder="Gen_molecular_pairs"):
    center_1, radius_1 = find_enclosing_sphere(coordinates_1)  # Find the center and radius of molecule 1
    center_2, radius_2 = find_enclosing_sphere(coordinates_2)  # Find the center and radius of molecule 2
    radius_2 -= radius_2_offset                                # Adjust radius of molecule 2 if needed

    # Generate random points on the surface of an ellipsoid defined by molecule 1
    random_coordinates_on_surface = generate_random_points(coordinates_1, radius_2, num_points)

    translated_molecules = []
    # Translate molecule 2 to various positions based on the random points
    for i in range(num_points):
        if i >= len(random_coordinates_on_surface):
            break
        translation_vector = random_coordinates_on_surface[i] - center_2  # Vector to translate molecule 2
        translated_coordinates = coordinates_2 + translation_vector       # Translate molecule 2
        translated_molecules.append(translated_coordinates)

    # Calculate the distance between molecule 1 and each translated molecule 2
    distances = [calculate_distance(coordinates_1, translated_coord) for translated_coord in translated_molecules]
    # Sort the molecules based on the distance
    sorted_distances_and_coordinates = sorted(zip(distances, translated_molecules), key=lambda x: x[0])

    merged_file_count = 0
    # Merge the molecules and write them to files
    for i, (distance, translated_coordinates) in enumerate(sorted_distances_and_coordinates):
        if merged_file_count >= 10:
            break
        if i % 2 == 0 or distance < 1.5:
            continue                    # Skip close molecules
        merged_file_count += 1
        merged_symbols = symbols_1 + symbols_2
        merged_coordinates = np.concatenate((coordinates_1, translated_coordinates))  # Combine the coordinates
        output_file_name = f"{file_name}_{merged_file_count}.xyz"
        output_file_path = os.path.join(output_folder, output_file_name)
        write_xyz_file(output_file_path, merged_symbols, merged_coordinates, distance)  # Write to file

    # If less than 10 files were created, try to add more
    if merged_file_count < 10:
        for i, (distance, translated_coordinates) in enumerate(sorted_distances_and_coordinates):
            if merged_file_count >= 10:
                break
            if i % 2 == 0 and distance >= 1.5:
                merged_file_count += 1
                merged_symbols = symbols_1 + symbols_2
                merged_coordinates = np.concatenate((coordinates_1, translated_coordinates))
                output_file_name = f"{file_name}_{merged_file_count}.xyz"
                output_file_path = os.path.join(output_folder, output_file_name)
                write_xyz_file(output_file_path, merged_symbols, merged_coordinates, distance)

    return merged_file_count
