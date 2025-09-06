from tkinter import font
import matplotlib.pyplot as plt
from ase.db import connect
import numpy as np
from sklearn.linear_model import LinearRegression

'''
This is to measure the number of atoms and hydrogen atoms in the database
'''

# List of database file paths and their display names
db_files = [
    ("MChem_DGMs/analysis/Databases/geom_drugs.db", "DRUGS"),
    ("MChem_DGMs/analysis/Databases/Gschnet_drugs_filtered.db", "Gschnet"),
    ("MChem_DGMs/analysis/Databases/EDM_drugs_filtered.db", "EDM"),
    ("MChem_DGMs/analysis/Databases/GeoLDM_drugs_filtered.db", "GeoLDM"),
    ("MChem_DGMs/analysis/Databases/JODO_drugs_filtered.db", "JODO")
]

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
