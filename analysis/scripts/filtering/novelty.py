from ase.db import connect
from ase.io import read, write
import json 


db = connect('/root/MChem_DGMs/analysis/Databases/JODO_qm9_cleaned_sorted_filtered.db')
smiles = json.load(open('/root/MChem_DGMs/analysis/Databases/qm9_smiles.json', 'r'))


non_novel_molecules = []
for row in db.select():
    if row.SMILES in smiles:
        non_novel_molecules.append(row)


print(f"Number of non-novel molecules: {len(non_novel_molecules)}")

