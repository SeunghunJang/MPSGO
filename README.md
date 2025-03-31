## Python framework for molecular pair structure generation and optimization (MPSGO)

Since chemical reactions essentially originate from molecular interactions, understanding molecular interactions beyond the properties of a single molecule is critical for the comprehensive investigation of real-world chemical experiments. However, quantum chemical calculation methods for molecule-molecule interaction need to consider extensive spatial configurations to find the optimal configuration. For this reason, a publicly available large database of molecular interaction results has not yet been established. Here, we devised an efficient and systematic framework for generating initial configurations of molecular interaction systems. Based on our framework, we constructed molecular interaction databases containing 49,620 individual molecules and 247,741 molecular interaction systems for chromophore-solvent, solute-solvent, and drug-drug molecular interactions. Our databases can be used for theoretical analysis and data-driven modeling of chemical reaction results, such as product property prediction and machine learning force fields. The molecular interaction databases will accelerate data-driven study in various research fields of chemical science.

## Run
- exec_smiles2xyz.py:         A Python script to convert SMILES strings in an Excel file into XYZ files for individual molecules.
- exec_xyz2molecular_pair.py: A python script to generate molecular pair structures using xyz files for individual molecules.
- exec_dft.py:                A python script for performing bulk quantum chemical calculations on the generated molecular pair structures.
- exec_extract_rlxed_xyz.py:  A python script to extract optimized molecular structure information from quantum chemistry calculation results and to generate xyz files.
- exec_extract_property.py:   A python script to extract quantum chemical property values ​​for optimized molecular structures from quantum chemical calculation results.

## Datasets
Any Excel file in <code style=" color: red">.xlsx</code> format containing SMILES string information is available for use, but the column headers for SMILES string information must be "SMILES_1" and "SMILES_2".
