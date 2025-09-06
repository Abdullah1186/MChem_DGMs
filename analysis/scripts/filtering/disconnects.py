
from ase.db import connect
from ase.io import read, write


db = connect('/root/MChem_DGMs/analysis/Databases/JODO_drugs_cleaned.db')


c=0
disconnected_molecules = []
for row in db.select():
    if '.' in row.orig_smiles: 
        c=c+1
        disconnected_molecules.append(row)

print(c)

