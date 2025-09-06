from ase.db import connect

# List of input database file paths
db_files = [
    "/home/chem/msummc/2192510_training_and_generation/drugs_subs/models/custom_data/drugs_run_30/generated_molecules/1.db"

    
]

# Path for the output concatenated database
output_db_path = "/home/chem/msummc/2192510_training_and_generation/drugs_subs/models/custom_data/drugs_run_30/generated_molecules/Gschnet_drugs_raw.db"

unique_entries = set()

# Create the output database and add unique entries
with connect(output_db_path, use_lock_file=False) as output_db:
    for db_file in db_files:
        with connect(db_file) as db:
            for row in db.select():
                # Define the unique property (e.g., chemical formula)
                if len(row.toatoms()) == 0:
                    continue
                
                unique_id = row.toatoms().get_chemical_formula()
                
                
                
                # Check if the unique ID has already been added
                if unique_id not in unique_entries:
                    output_db.write(row.toatoms(), key_value_pairs=row.key_value_pairs)
                    unique_entries.add(unique_id)