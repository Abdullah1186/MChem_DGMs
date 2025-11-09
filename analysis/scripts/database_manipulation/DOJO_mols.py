import pickle
from rdkit import Chem
from ase.db import connect
from ase import Atoms


def rdkit_mol_to_ase_atoms(rdkit_mol: Chem.Mol) -> Atoms:
    """Convert an RDKit molecule to an ASE Atoms object.
    
    Args:
        rdkit_mol: RDKit molecule object.
        
    Returns:
        ASE Atoms object.
    """
    ase_atoms = Atoms(
        numbers=[
            atom.GetAtomicNum() for atom in rdkit_mol.GetAtoms()
        ], 
        positions=rdkit_mol.GetConformer().GetPositions()
    )
    return ase_atoms


mol_list = pickle.load(open('/root/ancestral_ckpt_30_42.pkl', 'rb'))
db = connect('/root/JODO_qm9_raw.db')

with db:
    for mol in mol_list:
        atoms = rdkit_mol_to_ase_atoms(mol)
        smi = Chem.MolToSmiles(mol)
        db.write(atoms, key_value_pairs={'orig_smiles': smi})












