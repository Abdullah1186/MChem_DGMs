import os
from ase.io import read
from ase.db import connect

# Path to the folder containing .txt files
txt_folder = "/home/chem/msummc/EDM/outputs/edm_qm9_resume/eval/analyzed_molecules"

# Path to the output ASE database
output_db = "/home/chem/msummc/EDM/EDM.db"

# Initialize the ASE database
with connect(output_db, use_lock_file=False) as db:
    # Iterate through all files in the directory
    for file_name in os.listdir(txt_folder):
        if file_name.endswith(".txt"):
            txt_file = os.path.join(txt_folder, file_name)
            print(f"Processing: {txt_file}")
            try:
                # Read the .txt file into an ASE Atoms object
                atoms = read(txt_file, format="xyz")  # Use "xyz" format for XYZ-like data in .txt files

                # Check if valid Atoms object
                if not atoms or len(atoms) == 0:
                    print(f"Skipped {file_name}: No atoms found.")
                    continue

                # Get the chemical formula
                formula = atoms.get_chemical_formula()

                # Add to ASE database with metadata
                db.write(atoms, data={"chemical_formula": formula})
                print(f"Added {formula} from {file_name}.")
            except Exception as e:
                print(f"Failed to process {file_name}: {e}")

# Verify database
with connect(output_db) as db:
    count = len(list(db.select()))
    print(f"Database contains {count} entries.")