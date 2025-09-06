

# Write generated database to .mol_dict file
#python /root/MChem_DGMs/analysis/scripts/filtering/sort_db.py /root/MChem_DGMs/analysis/Databases/JODO_drugs_cleaned.db --write_mol_dict

# Filter .mol_dict file
python /root/MChem_DGMs/analysis/scripts/filtering/filter_generated.py  /root/MChem_DGMs/analysis/Databases/JODO_drugs_cleaned_sorted.mol_dict --filters 'valence' 'unique' --store 'valid' --print_file >> results.txt