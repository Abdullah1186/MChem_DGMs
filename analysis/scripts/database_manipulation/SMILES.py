import sqlite3
import pandas as pd
from ase.db import connect
import os
from rdkit import Chem
from IPython.display import display
import subprocess
import openbabel
import pybel
from rdkit.Chem import Draw
from tqdm import tqdm
import json 

# Define the base directory for databases
repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../.."))
base_dir = os.path.join(repo_root, "analysis", "Databases")


data_name = input("Enter the database name (with .db extension): ")

# Connect to the databas
db = connect(f"{base_dir}/{data_name}")

# Create a directory for structures if it doesn't exist
output_dir = f"{base_dir}/structures"
os.makedirs(output_dir, exist_ok=True)

dic={}
# Loop through all rows (molecules) in the database
for row in tqdm(db.select(), desc='loading'):  # Select all rows
    # Get the Atoms object from the current row
    molecule = row.toatoms()
  
    # Write the molecule to an .xyz file in the structures directory
    xyz_filename = os.path.join(output_dir, f'molecule_{row.id}.xyz')  # Save in 'structures'
    molecule.write(xyz_filename)

    # let pybel read file 
    mol = next(pybel.readfile("xyz", xyz_filename))
  
    # Convert to SMILES
    smiles = mol.write("smi").strip()
    smiles = smiles.split('\t')[0].strip()

    dic[f'molecule_{row.id}']=smiles
    



c=1

for i,row in tqdm(enumerate(db.select()),desc='updating'):
        db.update(id=row.id, SMILES= dic[f'molecule_{c}'])
        c=c+1




smiles_list = []
for row in db.select():
    smiles = row.get("SMILES")

    
    smiles_list.append(smiles)

new_data_name = data_name.split('.db')[0] 
f = open(f"{base_dir}/{new_data_name}_smiles.json", "w", encoding = "cp1252") # add your json file name
json.dump(smiles_list, f)
f.close()

print(f"SMILES have been added to the database and saved to {base_dir}/{new_data_name}_smiles.json")
