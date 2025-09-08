from ase.db import connect
from ase.io import read, write
import json 


db = connect('/root/MChem_DGMs/analysis/Databases/JODO_drugs_filtered.db')
smiles = json.load(open('/root/MChem_DGMs/analysis/Databases/geom_drugs_smiles.json', 'r'))


novel_molecules = []
for row in db.select():
    if row.SMILES in smiles:
        novel_molecules.append(row)


print(f"Number of novel molecules: {len(novel_molecules)}")

