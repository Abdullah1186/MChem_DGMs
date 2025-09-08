from tkinter import font
import matplotlib.pyplot as plt
from ase.db import connect
import numpy as np
import os 
from sklearn.linear_model import LinearRegression



'''
This is to measure the number of atoms and hydrogen atoms in the database
'''



script_dir = os.path.dirname(os.path.abspath(__file__))

# Repo root (go up 3 levels from SMILES.py in your case)
repo_root = os.path.abspath(os.path.join(script_dir, "../../.."))

# Databases directory
base_dir = os.path.join(repo_root, "analysis", "Databases")

choice1 = input(""" Choose database:
                1. QM9
                2. DRUGS
                3. OE62
                """)

choice2 = input(""" Use filtered database? (y/n):
                 """)

if choice1 == "1":
    if choice2 == "n":
        db_files = [
            (f"{base_dir}/Databases/qm9_smiles.db", "QM9"),
            (f"{base_dir}/Databases/Gschnet_qm9.db", "Gschnet"),
            (f"{base_dir}/Databases/EDM_qm9.db", "EDM"),
            (f"{base_dir}/Databases/GeoLDM_qm9.db", "GeoLDM")
        ]
    else:
        db_files = [
            (f"{base_dir}/qm9_smiles.db", "QM9"),
            (f"{base_dir}/Gschnet_qm9_filtered.db", "Gschnet_filtered"),
            (f"{base_dir}/EDM_qm9_filtered.db", "EDM_filtered"),
            (f"{base_dir}/GeoLDM_qm9_filtered.db", "GeoLDM_filtered")
        ]
elif choice1 == "2":
    if choice2 == "n":
        db_files = [
            (f"{base_dir}/geom_drugs.db", "DRUGS"),
            (f"{base_dir}/Gschnet_drugs_raw.db", "Gschnet"),
            (f"{base_dir}/EDM_drugs_raw.db", "EDM"),
            (f"{base_dir}/GeoLDM_drugs_raw.db", "GeoLDM"),
            (f"{base_dir}/JODO_drugs_raw.db", "JODO")
        ]
    else:
        db_files = [
            (f"{base_dir}/geom_drugs.db", "DRUGS"),
            (f"{base_dir}/Gschnet_drugs_filtered.db", "Gschnet_filtered"),
            (f"{base_dir}/EDM_drugs_filtered.db", "EDM_filtered"),
            (f"{base_dir}/GeoLDM_drugs_filtered.db", "GeoLDM_filtered"),
            (f"{base_dir}/JODO_drugs_filtered.db", "JODO_filtered")
        ]


elif choice1 == "3":
    if choice2 == "n":
        db_files = [
            (f"{base_dir}/OE62_full.db", "OE62"),
            (f"{base_dir}/Gschnet_oe62_raw.db", "Gschnet"),
            (f"{base_dir}/EDM_oe62_raw.db", "EDM")
        ]
    else:
        db_files = [
            (f"{base_dir}/OE62_full.db", "OE62"),
            (f"{base_dir}/Gschnet_oe62_filtered.db", "Gschnet_filtered"),
            (f"{base_dir}/EDM_oe62_filtered.db", "EDM_filtered")
        ]

else:
    quit()

fig, axes = plt.subplots(1, 5, figsize=(12, 3), sharey=True)

colors = ["black", "r", "b", "g","m"]


fig, axes = plt.subplots(1, 5, figsize=(20, 5), sharey=True, constrained_layout=True)

for ax, (db_path, title), color in zip(axes, db_files, colors):
    db = connect(db_path)
    num_atoms = []
    num_H_atoms = []

    for row in db.select():
        atoms = row.toatoms()
        symbols = atoms.get_chemical_symbols()
        num_atoms.append(len(symbols))
        num_H_atoms.append(symbols.count('H'))

    # Convert to numpy arrays
    X = np.array(num_atoms).reshape(-1, 1)
    y = np.array(num_H_atoms)

    # Fit linear regression
    model = LinearRegression()
    model.fit(X, y)
    slope = model.coef_[0]

    # Sorted line
    X_sorted = np.sort(X, axis=0)
    y_pred = model.predict(X_sorted)

    # Plot
    ax.scatter(X, y, s=10, alpha=0.5, c='black')
    ax.plot(X_sorted, y_pred, color=color, linewidth=2.5, label=f"Slope = {slope:.2f}")
    ax.set_title(title, fontsize=30)
    ax.set_xlabel("Atoms", fontsize=30)
    if ax == axes[0]:
        ax.set_ylabel("H atoms", fontsize=30)
    ax.legend(fontsize=20)
    ax.tick_params(axis='x', labelsize=20, )
    ax.tick_params(axis='y', labelsize=20)

plt.tight_layout()
plt.show()
