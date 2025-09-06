from ase.db import connect
from ase.io import read, write
from io import StringIO
from openbabel import pybel
import numpy as np

''' This script filters disconnected molecules in a database.
It retains the largest fragment and adds a SMILES representation of it.'''

db = connect('/root/MChem_DGMs/analysis/Databases/JODO_drugs_raw.db')
filtered_db_path = '/root/MChem_DGMs/analysis/Databases/JODO_drugs_cleaned.db'

def filter_disconnected(atoms):
    with StringIO() as fa:
        write(fa, atoms, format='xyz')
        mol = pybel.readstring('xyz', fa.getvalue())
        obmol = mol.OBMol
        frags = obmol.Separate()
        if len(frags) == 1:
            return None, None
        else:
            frag_na = [frag.NumAtoms() for frag in frags]
            largest_idx = np.argmax(frag_na)
            largest_mol = pybel.Molecule(frags[largest_idx])
            fb = StringIO(largest_mol.write('xyz'))
            return read(fb, format='xyz'), largest_mol.write('can').strip()

with connect(filtered_db_path) as db_filtered:
    for i, row in enumerate(db.select()):
        atoms = row.toatoms()
        filtered_atoms, smi = filter_disconnected(atoms)
        if filtered_atoms is None:
            db_filtered.write(atoms, key_value_pairs=row.key_value_pairs, data=row.data)
        else: 
            kvp = row.key_value_pairs
            kvp['orig_smiles'] = smi
            db_filtered.write(filtered_atoms, key_value_pairs=kvp, data=row.data)