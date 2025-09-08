import itertools
from tkinter import font
import numpy as np
import matplotlib.pyplot as plt
from ase.db import connect
from ase.neighborlist import natural_cutoffs, neighbor_list
from scipy.stats import gaussian_kde
from elemental_composition_analysis import sample_generated
import matplotlib as mpl
import os 


script_dir = os.path.dirname(os.path.abspath(__file__))

# Repo root (go up 3 levels from SMILES.py in your case)
repo_root = os.path.abspath(os.path.join(script_dir, "../../.."))

# Databases directory
base_dir = os.path.join(repo_root, "analysis", "Databases")

size=20
params = {'legend.fontsize': size,
          'figure.figsize': (7,5),
          'lines.linewidth': 5,
          'axes.labelsize': size,
          'axes.titlesize': size,
          'xtick.labelsize': size,
          'ytick.labelsize': size,
          'axes.titlepad': 10,
         'axes.linewidth':1.5,
         'axes.prop_cycle': mpl.cycler(color=["k","r", "b", "g"])
                                       }
plt.rcParams.update(params)

def set_label(database):
    label_db=f"{database}"
    if database=="qm9_smiles.db":
            label_db="QM9"
    elif "GeoLDM" in database:
            label_db="GeoLDM"
    elif "EDM" in database:
            label_db="EDM"#"Generated"
    elif "geom" in database:
            label_db="GEOM-DRUGS"#"Generated"
    # elif "qm9" in database:
    #         label_db="QM9"
    elif "OE62_full" in database:
            label_db="OE62" 
            
    elif "Gschnet" in database:
            label_db="Gschnet"


    
  #  elif "filtered"
    return label_db

def get_interatomic_distances(database, atom_pair, sample = None):
    """Returns a list of the distances between atoms matching the given pair
    for each molecule in the given database.
    """

    distances_list = []

    db = connect(f"{base_dir}/{database}")

    # Convert rows to atoms and get the chemical symbols of atoms in the 
    # molecule.
    for row in db.select():
        if sample is not None:
            if row.id - 1 not in sample:
                continue
        atoms = row.toatoms()
        symbols = atoms.get_chemical_symbols()
        indices = []

        # Get the index of each atom paired with its symbol if that symbol
        # matches the atom pair given.
        for i, symbol in enumerate(symbols):
            if symbol in atom_pair:
                indices.append((i, symbol))

        molecule_distances = []

        # Here, i represents the index of the index-symbol pair in the indices
        # list. j is the index-symbol pair. If the two atoms in the pair given
        # are the same, the distances between the atom represented by j and all
        # other matching atoms is obtained. If not, a list containing all
        # indices of atoms with elements not matching that of the atom
        # represented by i is created. n represents these indices. The distance
        # between atom i and these atoms are then obtained.
        for i, j in enumerate(indices[:-1]):
            if atom_pair[0] == atom_pair[1]:
                distances = atoms.get_distances(j[0], [j[0] for j in indices[i+1:]])
            else:
                opposite_indices = [j[0] for j in indices[i+1:] if j[1] != indices[i][1]]
                if len(opposite_indices) == 0:
                    continue
                distances = atoms.get_distances(j[0], [n for n in opposite_indices])
            for distance in distances:
                molecule_distances.append(distance)

        distances_list.append(molecule_distances)

    return distances_list


def plot_distances(databases, atom_pair, samples = None, bonded = False):
    """Plots a histogram showing the distribution of distances between atoms
    matching the given pair for molecules found in the given database.
    """
    for i, database in enumerate(databases):
        if bonded == False:
            if i == 0 or samples == None:
                distances_list = get_interatomic_distances(database[0], atom_pair)
            else:
                distances_list = get_interatomic_distances(database[0], atom_pair, sample = samples[i - 1])
            all_distances = np.array(list(itertools.chain.from_iterable(distances_list)))
        else:
            molecules = get_neighbours(database)
            all_distances = []
            with connect(f"{base_dir}/{database[0]}") as db:
                for idx, row in enumerate(db.select()):
                    atoms = row.toatoms()
                    atom1s = molecules[idx][0]
                    atom2s = molecules[idx][1]
                    for atom_idx in range(0, len(atom1s)):
                        if {atoms[atom1s[atom_idx]].symbol, atoms[atom2s[atom_idx]].symbol} == set(atom_pair):
                            if i == 0 or samples == None:
                                all_distances.append(molecules[idx][2][atom_idx])
                            else:
                                if idx in samples[i - 1]:
                                    all_distances.append(molecules[idx][2][atom_idx])
            all_distances = np.array(all_distances)


        density = gaussian_kde(all_distances, bw_method = 0.1)
        x = np.linspace(0, all_distances.max(), 1000)
        plt.plot(x, density(x), label = set_label(database[0]))

    plt.xlabel("Distance (Å)",fontsize=40)
    plt.ylabel("Density",fontsize=40)
    plt.xlim(1,1.7)
    # if len(databases) == 1:
    #     plt.title(f"KDE Showing Distribution of {atom_pair[0]}-{atom_pair[1]} Distances in {databases[0]}")
    # elif len(databases) == 2:
    #     plt.title(f"KDE Showing Distribution of {atom_pair[0]}-{atom_pair[1]} Distances in QM9 Databases")
    # elif len(databases) == 3:
    #     plt.title(f"KDE Showing Distribution of {atom_pair[0]}-{atom_pair[1]} Distances in OE62+THz Databases")
    # elif len(databases) == 5:
    #     plt.title(f"KDE Showing Distribution of {atom_pair[0]}-{atom_pair[1]} Distances in OE62 Databases")
    # else:
    #     plt.title(f"KDE Showing Distribution of {atom_pair[0]}-{atom_pair[1]} Distances in Thiols Databases")
    plt.margins(x=0)
    plt.legend(fontsize=40)
    plt.tight_layout()
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.show()


def get_neighbours(database):
    """Return a list of neighbour lists with the two atoms and the distance
    between them for a single database.
    """
    neighbours = []
    with connect(f"{base_dir}/{database[0]}") as db:
        for row in db.select():
            atoms = row.toatoms()
            cutoffs = natural_cutoffs(atoms, mult = 1.1)
            i, j, d = neighbor_list("ijd", atoms, cutoffs)
            molecule_neighbours = (i, j, d)
            neighbours.append(molecule_neighbours)
    
    return neighbours


def show_extreme_neighbours(databases):
    """Output the nearest and furthest neighbours in a database."""

    for database in databases:
        max_distances = []
        min_distances = []
        i = 0
        molecules = get_neighbours(database)
        for molecule in molecules:
            distances = molecule[2]
            try:
                max_distances.append(np.max(distances))
                min_distances.append(np.min(distances))
            except ValueError:
                print(i)
            finally:
                i += 1

        print(f"Furthest neighbours in {database[0]} are {max(max_distances)} Å apart.")
        print(f"Nearest neighbours in {database[0]} are {min(min_distances)} Å apart.")


def main():
    """The main program, where the user enters inputs to choose which 
    functions to run
    """

    while True:
        # Creating argument sets to pass into the functions
        # to prevent input errors, based on a choice.
        choice1 = input(
            """Enter the number corresponding to the database you would like to use
            (enter anything else to end):
            1. QM9
            2. DRUGS
            3. OE62
            4. DATABASES
            """
            )

        filtered = input("Use filtered data? (y/n)\n")


        if choice1 == "1": 
            if filtered == 'n':## raw 
                databases = (
                ("qm9_smiles.db", "qm9_ecomp.json"),
                ("Gschnet_qm9.db", "Gschnet_qm9_ecomp.json"),
                
                ("EDM_qm9.db","EDM_qm9_ecomp.json"),
                ("GeoLDM_qm9.db", "GeoLDM_qm9_ecomp.json")
                )
            else:
                databases = (
               ("qm9_smiles.db", "qm9_ecomp.json"),
               ("Gschnet_qm9_filtered.db", "Gschnet_qm9_filtered_ecomp.json"),
               ("EDM_qm9_filtered.db","EDM_qm9_filtered_ecomp.json"),
               ("GeoLDM_qm9_filtered.db","GeoLDM_qm9_filtered_ecomp.json")
               )##filtered

        elif choice1 =="2":
            if filtered == 'y':##raw6
               databases = (
                    ("geom_drugs.db","geom_drugs_ecomp.json"),
                    ("Gschnet_drugs_filtered.db","Gschnet_drugs_filtered_ecomp.json"),
                
                    ("EDM_drugs_filtered.db","EDM_drugs_filtered_ecomp.json"),
                    ("GeoLDM_drugs_filtered.db","GeoLDM_drugs_filtered_ecomp.json")
                )
            else:
                databases = (
                    ("geom_drugs.db","geom_drugs_ecomp.json"),
                    ("Gschnet_drugs_raw.db","Gschnet_drugs_raw_ecomp.json"),

                    ("EDM_drugs_raw.db","EDM_drugs_raw_ecomp.json"),
                    ("GeoLDM_drugs_raw.db","GeoLDM_drugs_raw_ecomp.json"),
                    ("JODO_drugs_raw.db","JODO_drugs_raw_ecomp.json")
                )


        elif choice1 == "3":
            if filtered == 'n':
                databases = (
                ("OE62_full.db", "OE62_full_ecomp.json"),
                    ("Gschnet_oe62_raw.db", "Gschnet_oe62_raw_ecomp.json"),
                    ("EDM_oe62_raw.db", "EDM_oe62_raw_ecomp.json")

                    )
                
            else:
                databases = (
                    ("OE62_full.db", "OE62_full_ecomp.json"),
                    ("Gschnet_oe62_filtered.db", "Gschnet_oe62_filtered_ecomp.json"),
                    ("EDM_oe62_filtered.db", "EDM_oe62_filtered_ecomp.json")
                )



        elif choice1 == "4":
            databases = (
               ("OE62_full.db", "OE62_full_ecomp.json"),
               ("qm9_smiles.db", "qm9_ecomp.json"),
                ("geom_drugs.db","geom_drugs_ecomp.json")
            )
            
        

        

        else:
            quit()

        choice2 = input(
            """What would you like to do? Enter the corresponding number, or 
            anything else to quit.
            1. Plot atoms distances.
            2. Get nearest and furthest neighbours.
            """
        )

        if choice2 == "1":
            #try:
            atom_pair = input(
                """Enter the pair of atoms you would like to measure
                the distance between, in the format A-B, where A and B
                are elements (they can be the same)
                """
                ).split("-")
            bonded = False
            bonded_only = input("Consider only bonded atoms? (y/n)\n")
            if bonded_only == "y":
                bonded = True
            choice3 = input("Use sampled generated databases? (y/n)\n")
            if choice3 == "y":
                select_data = input(
                    """What would you like to sample according to?
                    1. Molecular Weight
                    2. Number of Atoms
                    """
                    )
                if select_data == "1":
                    data = "Molecular Weights"
                    bin_width = 10
                elif select_data == "2":
                    data = "Number of Atoms"
                    bin_width = 1
                else:
                    continue
                generated_samples = sample_generated(databases, sample_data = data, bin_width = bin_width)
                plot_distances(databases, atom_pair, generated_samples, bonded = bonded)
            else:
                plot_distances(databases, atom_pair, bonded = bonded)
            #except ValueError:
                #print("No molecules satisfied the conditions.")
        elif choice2 == "2":
            show_extreme_neighbours(databases)

if __name__ == "__main__":
    main()

