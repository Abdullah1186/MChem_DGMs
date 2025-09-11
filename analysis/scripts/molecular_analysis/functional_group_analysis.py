"""This will be used to generate plots showing the distribution of
functional groups throughout the databases, and other miscellaneous
data.
"""
import pickle
import os
import sys
import json
import itertools
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("TkAgg")
from rdkit import Chem
from rdkit.Chem import RDConfig
from rdkit.Chem import Lipinski
from rdkit import RDLogger
from ase.db import connect
from elemental_composition_analysis import get_valid, sample_generated
import matplotlib as mpl


sys.path.append(os.path.join(RDConfig.RDContribDir, 'IFG'))
from ifg import identify_functional_groups

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
         'axes.prop_cycle': mpl.cycler(color=["k","r", "b","g","y"])
                                       }
plt.rcParams.update(params)


def set_label(database):
    label_db=f"{database}"

    if "qm9_smiles.db" in database:
        label_db = "QM9"


    elif 'geom' in database:
        label_db="GEOM-DRUGS"

    elif "OE62_full" in database:
        label_db="OE62"

    elif "Gschnet" in database:
        label_db="Gschnet"

    elif "JODO" in database:
        label_db="JODO"

    elif "EDM" in database:
        label_db="EDM"

    elif "GeoLDM" in database:
        label_db="GeoLDM"

    if "filtered" in database:
        label_db += "_filtered"
    else: 
        label_db += ""

    return label_db

def generate_smiles(database, smiles_file):
    """Defining a function to generate smiles data for a chosen database."""

    db = connect(f"{base_dir}/{database}")
    smiles_list = []

    for row in db.select():
        smiles = row.get("SMILES") or row.get("smiles")
        if smiles == "":
            smiles = "INVALID"
        
        smiles_list.append(smiles)

    f = open(f"{base_dir}/{smiles_file}", "w", encoding = "cp1252")
    json.dump(smiles_list, f)
    f.close()


def find_substructure(database, smiles_file, substructure,sample):
    """This will search the smiles file for molecules containing the
    desired substructure.
    """

    f = open(f"{base_dir}/{smiles_file}", "r", encoding = "cp1252")
    smiles_list = json.load(f)
    f.close()
    patt = Chem.MolFromSmarts(substructure)
    matches = []
    smiles_matches = []
    failed = []
    
## new code ("db","smile")
    # if database in ["qm9_filtered.db", "OE62_full.db","geom_drugs.db","qm9_smiles.db"]:
    #     new_smiles_list = smiles_list

    if sample is not None:
        new_smiles_list= [samp for i, samp in enumerate(smiles_list) if i in sample]

    else:
        new_smiles_list = smiles_list
       
        #if i!=0 sample[i-1]
## newcode 
    for smiles in new_smiles_list:
        try:
            molecule = Chem.MolFromSmiles(smiles)
            if molecule.HasSubstructMatch(patt):
                matches.append(molecule)
                smiles_matches.append(smiles)
        except AttributeError:
            failed.append(smiles)

    db = connect(f"{base_dir}/{database}")

    if database in ["qm9_filtered.db", "OE62_full.db","geom_drugs.db","qm9_smiles.db"]: ## data bases that are training data
        #molecules_count = db.count()
        molecules_count = len(new_smiles_list)
        print(f"{len(matches)} matches found in {database} out of {molecules_count} molecules ({(len(matches)/molecules_count)*100:.2f}%), of which {len(failed)} RDKit failed to convert (considering only valid molecules, {(len(matches)/(molecules_count-len(failed)))*100:.2f}%)")

    else: 
        molecules_count = len(new_smiles_list)
        print(f"{len(matches)} matches found in {database} out of {molecules_count} molecules ({(len(matches)/molecules_count)*100:.2f}%), of which {len(failed)} RDKit failed to convert (considering only valid molecules, {(len(matches)/(molecules_count-len(failed)))*100:.2f}%)")


    return matches, smiles_matches


def get_appearances(database, smiles_file, substructure):
    """Gets the number of molecules where the specified substructure
    appears.
    """

    matches, _ = find_substructure( database, smiles_file, substructure)
    apperances = len(matches)
    return apperances


def get_instances(database, smiles_file, substructure,sample = None):
    """Returns a dictionary containing the number of times that each
    """

    matches, _ = find_substructure(database, smiles_file, substructure,sample=sample)
    instances = []
    patt = Chem.MolFromSmarts(substructure)

    for molecule in matches:
        count = len(molecule.GetSubstructMatches(patt))
        instances.append(count)

    frequencies_set = sorted(list(set(instances)))
    counts = [instances.count(times) for times in frequencies_set]
    instance_dict = {instances_in_molecule:count for instances_in_molecule, count in zip(frequencies_set, counts)}

    return instance_dict


def create_barplot(database, smiles_file, substructure, group_name):
    """A function which plots a barplot showing how many times a functional
    group appears in molecules across the database.
    """

    instance_dict = get_instances(database, smiles_file, substructure)
    instances = list(map(str, list(instance_dict.keys())))
    counts = np.array(list(instance_dict.values()), dtype = float)
    counts /= len(get_valid(smiles_file))

    plt.bar(instances, counts, edgecolor = "black")

    plt.xlabel("Instances of Substructure in Molecule")
    plt.ylabel("Proportion of Total Valid Molecules")
    plt.title(f"Barplot Showing Instances of {group_name} Per Molecule in {database}")
    # plt.show()


def create_barplot_all(args, substructure, group_name,sample = None):
    """A function which plots a barplot showing how many times a functional
    group appears in molecules across all databases.
    """
    ## new code to sample generated data based on number of atoms of training data 
    # list of sampled id's 
   #generated_samples = sample_generated(args, sample_data = "Number of Atoms", bin_width = 1, use_valid = False, split = False )
    
    
    

    ## 
    num_datasets = len(args)
    nos_of_instances = []
    proportions_list = []

    for idx,arg_set in enumerate(args):
        if idx ==0 or sample is None:
            instance_dict = get_instances(arg_set[0] , arg_set[1], substructure)

        else:
            instance_dict = get_instances(arg_set[0] , arg_set[1], substructure,sample=sample[idx-1])
        
        instances = list(map(str, list(instance_dict.keys())))
        #new code only taking instances that are in sampled 
        
    
        
        for inst in instances:
            if inst not in nos_of_instances:
                nos_of_instances.append(inst)

        counts = np.array(list(instance_dict.values()), dtype = float)
        if idx == 0 or sample is None:
            counts /= len(get_valid(arg_set[1]))
        else:
            counts /= len(sample[idx-1])
        proportions_list.append(counts)

    nos_of_instances.sort(key = int)

    x = list(range(1, len(nos_of_instances) + 1))
    bar_width = 0.1
    shifts = np.arange(-num_datasets / 2, num_datasets / 2) * bar_width + bar_width/2
    shifted_xs = [x + shifts[i] for i in range(num_datasets)]

    # Correcting the length of the proportions arrays if they do not match
    # the length of nos_of_instances, and therefore shifted_xs
    proportions_list_array = np.array(proportions_list, dtype=object)
    # np.save(f"/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/shifted_xs_{group_name}.npy",shifted_xs)
    # np.save(f"/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/proportion_list_{group_name}.npy",proportions_list_array)
    # np.save(f"/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/nos_of_instances_{group_name}.npy",nos_of_instances) 

    for i in range(num_datasets):
        if len(proportions_list[i]) < len(nos_of_instances):
            new_proportions = np.zeros(len(nos_of_instances), dtype = float)
            new_proportions[:len(proportions_list[i])] = proportions_list[i]
            proportions_list[i] = new_proportions

        plt.bar(
            shifted_xs[i],
            proportions_list[i],
            bar_width,
            label = set_label(args[i][0]),
            edgecolor = "black"
            )

    plt.xlabel("Instances of Substructure in Molecule", fontsize=40)
    plt.ylabel("Proportion of Total Molecules", fontsize=40)
    plt.title(f"Barplot Showing Instances of {group_name} Per Molecule ")
    plt.tick_params(axis='both', labelsize=30)
    plt.xticks(np.arange(1, len(x) + 1), labels = nos_of_instances)
    plt.margins(x = 0)
    plt.legend(fontsize=30)
    plt.show()
    plt.savefig("plot.png", dpi=300)



def get_aro_member_counts(mol, ring_info):
    """Gets the number of atoms that appear in aromatic rings in a molecule.
    """

    aro_ring_sizes = []
    rings_bonds = ring_info.BondRings()
    
    for ring in rings_bonds:
        is_aromatic = lambda bond_idx : mol.GetBondWithIdx(bond_idx).GetIsAromatic()
        if all(list(map(is_aromatic, ring))):
            aro_ring_sizes.append(len(ring))
    
    return aro_ring_sizes if len(aro_ring_sizes) > 0 else None


def get_misc_info(database, smiles_file):
    """Gets information such as the number of rings, or the degree of
    unsaturation (or index of hydrogen deficiency, IDH)
    """

    f = open(f"{base_dir}/{smiles_file}", "r", encoding = "cp1252")
    smiles_list = json.load(f)
    f.close()

    ring_nums = []
    ring_type_counts = {"aromatic": 0, "unsaturated aliphatic": 0, "saturated": 0}
    aro_rings_atom_counts = []
    ring_type_nums = {"aromatic": [], "unsaturated aliphatic": [], "saturated": []}
    ihd_list = []

    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            ring_info = mol.GetRingInfo()
            ring_nums.append(ring_info.NumRings())
            aro_member_count = get_aro_member_counts(mol, ring_info)
            if aro_member_count is not None:
                aro_rings_atom_counts.append(aro_member_count)

            num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "C")
            num_hydrogens = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms())
            num_nitrogens = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "N")
            num_halogens = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() in ["F", "Cl", "Br", "I"])

            degree_of_unsaturation = (2 * num_carbons + 2 - num_hydrogens + num_nitrogens - num_halogens) / 2
            ihd_list.append(degree_of_unsaturation)

            num_aromatic = Lipinski.NumAromaticRings(mol)
            num_aliphatic = Lipinski.NumAliphaticRings(mol)
            num_saturated = Lipinski.NumSaturatedRings(mol)
            num_unsaturated = num_aliphatic - num_saturated

            ring_type_nums["aromatic"].append(num_aromatic)
            if num_aromatic > 0:
                ring_type_counts["aromatic"] += 1
            ring_type_nums["unsaturated aliphatic"].append(num_unsaturated)
            if num_unsaturated > 0:
                ring_type_counts["unsaturated aliphatic"] += 1
            ring_type_nums["saturated"].append(num_saturated)
            if num_saturated > 0:
                ring_type_counts["saturated"] += 1

        except AttributeError:
            continue
    
    aro_rings_atom_counts = list(itertools.chain.from_iterable(aro_rings_atom_counts))

    try:
        avg_ring_num = sum(ring_nums)/len(ring_nums)
    except ZeroDivisionError:
        avg_ring_num = 0
    try:
        avg_ihd = sum(ihd_list)/len(ihd_list)
    except ZeroDivisionError:
        avg_ihd = 0

    print(f"Average number of rings per molecule in {database} is {avg_ring_num:.2f}")
    print(f"Average index of hydrogen deficiency (IHD) per molecule in {database} is {avg_ihd:.2f}")

    return ring_nums, ring_type_nums, aro_rings_atom_counts, ring_type_counts, avg_ihd


def plot_misc_data(args, rings, ring_type_nums, aro_atoms, ring_type_counts, ihd):
    """A function for plotting the distribution of rings across
    databases in a histogram and the average IHD of each database in a
    barplot.
    """

    # IHD barplot
    # num_datasets = len(args)
    # databases = [arg_set[0] for arg_set in args]
    # plt.bar(databases, ihd, edgecolor = "black", align = "edge", width = 0.3)
    # database_labels = ["\n" + set_label(database) if databases.index(database) % 2 == 0 else set_label(database) for database in databases]

    # plt.xlabel("Database")
    # plt.ylabel("Average IHD")
    # if len(args) == 8:
    #     plt.title("Barplot Showing Average IHD of molecules across Thiols databases")
    # elif len(args) == 5:
    #     plt.title("Barplot Showing Average IHD of molecules across OE62 databases")
    # elif len(args) == 3:
    #     plt.title("Barplot Showing Average IHD of molecules across OE62+THz databases")
    # else:
    #     plt.title("Barplot Showing Average IHD of molecules across qm9 databases")
    # plt.xticks(np.arange(num_datasets), labels = database_labels)
    # plt.margins(x = 0)
    # plt.show()

    #Total rings barplot
    ring_types = ("none", "aromatic", "unsaturated aliphatic", "saturated")
    # num_datasets = len(args)
    # bar_width = 0.8 / num_datasets
    # shifts = np.arange(-num_datasets / 2, num_datasets / 2) * bar_width + bar_width/2

    # plt.figure(figsize=(20, 20))

    # for ring_type in ring_types:
    #     for i, arg_set in enumerate(args):
    #         if ring_type == "none":
    #             rings_data = rings[i]
    #         else:
    #             rings_data = ring_type_nums[i][ring_type]
    #         bins = np.arange(max(rings_data) + 1) - 0.5
    #         shifted_bins = bins + shifts[i]
    #         plt.hist(
    #             rings_data,
    #             bins = shifted_bins,
    #             edgecolor = "black",
    #             width = bar_width,
    #             rwidth = bar_width,
    #             label = set_label(arg_set[0]),
    #             density = True
    #             )

    #     if ring_type == "none":
    #         ring_type = "all"

    #     plt.xlabel("Count")
    #     plt.ylabel("Density")
    #     if len(args) == 8:
    #         plt.title(f"Histogram Showing Distribution of {ring_type} Ring Numbers in Thiols databases")
    #     elif len(args) == 5:
    #         plt.title(f"Histogram Showing Distribution of {ring_type} Ring Numbers in OE62 databases")
    #     elif len(args) == 3:
    #         plt.title(f"Histogram Showing Distribution of {ring_type} Ring Numbers in OE62")
    #     else:
    #         plt.title(f"Histogram Showing Distribution of {ring_type} Ring Numbers in qm9 databases")
    #     plt.margins(x=0)
    #     if ring_type == "all":
    #         plt.xticks(np.arange(max([max(rings_data) for rings_data in rings]) + 1))
    #     else:
    #         plt.xticks(np.arange(max([max(rings_data[ring_type]) for rings_data in ring_type_nums]) + 1))
    #     plt.legend()
        #plt.show()

    # Aromatic rings barplot
    # num_atoms = []

    # for database in aro_atoms:
    #     num_atoms += list(set(database))

    # num_atoms = sorted(list(set(num_atoms)))
    # xs = np.arange(len(num_atoms))
    # x_labels = [str(i) for i in num_atoms]
    # shifts = np.arange(-num_datasets / 2, num_datasets / 2) * bar_width + bar_width/2

    # for i, arg_set in enumerate(args):
    #     heights = []
    #     for num in num_atoms:
    #         heights.append(aro_atoms[i].count(num))
    #     heights = (np.array(heights)/len(aro_atoms[i])) * 100
    #     shifted_xs = xs + shifts[i]

    #     plt.bar(
    #         shifted_xs,
    #         heights,
    #         width = bar_width,
    #         edgecolor = "black",
    #         label = set_label(arg_set[0])
    #     )
    
    # plt.xticks(np.arange(len(num_atoms)), labels = x_labels)
    # plt.xlabel("Atoms in Ring")
    # plt.ylabel("Density")
    # if len(args) == 8:
    #     plt.title("Barplot Showing Number of Atoms in Aromatic Rings in Thiols databases")
    # elif len(args) == 5:
    #     plt.title("Barplot Showing Number of Atoms in Aromatic Rings in OE62 databases")
    # elif len(args) == 3:
    #     plt.title("Barplot Showing Number of Atoms in Aromatic Rings in OE62+THz databases")
    # else:
    #     plt.title("Barplot Showing Number of Atoms in Aromatic Rings in qm9 databases")
    # plt.legend()
    # plt.margins(x=0)
    # plt.show()

    # Ring Types barplot
    # xs = np.arange(len(ring_types))
    # x_labels = ring_types
    # shifts = np.arange(-num_datasets / 2, num_datasets / 2) * bar_width + bar_width/2

    # for i, arg_set in enumerate(args):
    #     heights = []
    #     for ring_type in ring_types:
    #         if ring_type == "none":
    #             heights.append((rings[i].count(0) / len(get_valid(arg_set[1]))) * 100)
    #         else:
    #             heights.append((ring_type_counts[i][ring_type] / len(get_valid(arg_set[1]))) * 100)
    #     shifted_xs = xs + shifts[i]

    #     plt.bar(
    #         shifted_xs,
    #         heights,
    #         width = bar_width,
    #         edgecolor = "black",
    #         label = set_label(arg_set[0])
    #     )

    # plt.xticks(np.arange(len(ring_types)), labels = x_labels)
    # plt.xlabel("Ring Type")
    # plt.ylabel("Density")
    # if len(args) == 8:
    #     plt.title("Barplot Showing Proportion of Molecules with Different Ring Types in Thiols databases")
    # elif len(args) == 5:
    #     plt.title("Barplot Showing Proportion of Molecules with Different Ring Types in OE62 databases")
    # elif len(args) == 3:
    #     plt.title("Barplot Showing Proportion of Molecules with Different Ring Types in OE62")
    # else:
    #     plt.title("Barplot Showing Proportion of Molecules with Different Ring Types in qm9 databases")
    # plt.legend()
    # plt.margins(x=0)
    # plt.show()

    #Combined types/counts (single database)
    colors = ["blue", "orange", "green", "red", "purple"]
    split = False
    # if "OE62" in args[0][0] and len(args) == 2:
    #     args[0][0] = "OE62"
    #     args[1][0] = "Gen."
    # elif "OE62" in args[0][0] and len(args) == 3:
    #     args[0][0] = "OE62"
    #     args[1][0] = "Raw Gen."
    #     args[2][0] = "Gen."
    # elif len(args) == 3:
    #     split_data = input("Split database into OE62 and THz (y/n)?\n")
    #     if split_data == "y":
    #         split = True
    #     if split == True:
    #         args[0][0] = "OE62"
    #         args.insert(1, ("THz", args[0][1]))
    #         args[2][0] = "Generated\n(transform)"
    #         args[3][0] = "Generated\n(no transform)"
    #     else:
    #         args[0][0] = "OE62+THz"
    #         args[1][0] = "Generated\n(transform)"
    #         args[2][0] = "Generated\n(no transform)"
    
    # elif len(args) == 2:
    #     args[0][0] = "QM9"
    #     args[1][0] = "Gen."


    # elif len(args) == 5:
    #     args[0][0] = "OE62"
    #     args[1][0] = "Raw Gen 1"
    #     args[2][0] = "Filtered Gen 1"
    #     args[3][0] = "Raw Gen 2"
    #     args[4][0] = "Filtered Gen 2"

    bar_width = 0.8 / len(args)
    shifts = np.arange(-len(args) / 2, len(args) / 2) * bar_width + bar_width/2

    use_sampled = input("Use sampled generated databases? (y/n)\n")
    if use_sampled == "y":
        select_data = input(
            """What would you like to sample according to?
            1. Molecular Weight
            2. Number of Atoms
            """
            )
        if select_data == "1":
            data = "Molecular Weights"
            bin_width = 10
        else:
            data = "Number of Atoms"
            bin_width = 1
        original_args = [arg_set for arg_set in args]
        # if len(args) == 4:
        #     args.insert(0, ["OE62+THz_natoms_60.db", "ot60_smiles.json"])
        #     del args[2]
        #     del args[1]
        generated_samples = sample_generated(args, sample_data = data, bin_width = bin_width, use_valid = True, split = split)
        args = original_args

    for i, arg_set in enumerate(args):
        database = arg_set[0]

        height_none = 0
        heights_aromatic = {}
        heights_unsaturated = {}
        heights_saturated = {}

        contains_gold = []
        valid = get_valid(arg_set[1])
        if valid[-1] == "INVALID":
            valid = valid[:-1]

        

        for ring_type in ring_types:
            if ring_type == "none":

                if use_sampled == "y" and i != 0:
                    rings_list = [rings for idx, rings in enumerate(rings[i]) if idx in generated_samples[i - 1]]
                    height_none = (rings_list.count(0)) / len(generated_samples[i - 1])
                else:
                    rings_list = rings[i]
                    height_none = (rings_list.count(0)) / len(get_valid(arg_set[1]))

            else:
                rings_data = ring_type_nums[i][ring_type]
                stack_count = ring_type_counts[i][ring_type]
                if use_sampled == "y":
                    rings_data = [data for idx, data in enumerate(rings_data) if idx in generated_samples[i - 1]]
                else: 
                    rings_data = ring_type_nums[i][ring_type]
                    
                

                if use_sampled == "y":
                    stack_density = stack_count / len(generated_samples[i - 1])
                else:
                    stack_density = stack_count / len(get_valid(arg_set[1])) 

                    
                unique_nums = sorted(list(set(rings_data)))

                for j in unique_nums:
                    if ring_type == "aromatic":
                        heights_aromatic[j] = (rings_data.count(j) / stack_count) * stack_density
                    elif ring_type == "unsaturated aliphatic":
                        heights_unsaturated[j] = (rings_data.count(j) / stack_count) * stack_density
                    elif ring_type == "saturated":
                        heights_saturated[j] = (rings_data.count(j) / stack_count) * stack_density
    
        xs = np.arange(len(ring_types))
        shifted_xs = xs + shifts[i]
        x_labels = ring_types
        aromatic_keys = list(heights_aromatic.keys())
        unsaturated_keys = list(heights_unsaturated.keys())
        saturated_keys = list(heights_saturated.keys())
        all_keys = set(aromatic_keys + unsaturated_keys + saturated_keys)

        for j in all_keys:
            if j not in aromatic_keys:
                heights_aromatic[j] = 0
            if j not in unsaturated_keys:
                heights_unsaturated[j] = 0
            if j not in saturated_keys:
                heights_saturated[j] = 0

        new_keys = ["1", "2", "3", "4+"]

        def group_4_or_greater(heights_dict):
            values = list(heights_dict.values())
            new_values = [values[i] for i in range(1, 4)]
            new_values.insert(3, sum(values[4:]))
            new_heights_dict = {key:value for key, value in zip(new_keys, new_values)}
            return new_heights_dict
        
        heights_aromatic = group_4_or_greater(heights_aromatic)
        heights_unsaturated = group_4_or_greater(heights_unsaturated)
        heights_saturated = group_4_or_greater(heights_saturated)

        heights_list = np.array([height_none, 0.0, 0.0, 0.0])
        heights_bottoms = np.array([0.0, 0.0, 0.0, 0.0])
        new_keys.append("extra")

        for j, key in enumerate(new_keys):
            color = colors[j]
            
            if j == 4:
                j = "4+"
            if i == 0:
                plt.bar(
                    shifted_xs,
                    heights_list,
                    bottom = heights_bottoms,
                    width = bar_width,
                    color = color,
                    edgecolor = "black",
                    label = f"{j} rings"
                )
            else:
                plt.bar(
                    shifted_xs,
                    heights_list,
                    bottom = heights_bottoms,
                    width = bar_width,
                    color = color,
                    edgecolor = "black"
                )
            heights_bottoms += heights_list
            try:
                heights_list = [0, heights_aromatic[key], heights_unsaturated[key], heights_saturated[key]]
            except KeyError:
                continue
        
        for x in shifted_xs:
            #plt.text(x, -0.22, set_label(f"{database}"), ha = "center", va = "center", rotation = 45, fontsize = 20)
            plt.text(x, -0.2, set_label(f"{database}"),ha="center", va="center", rotation=90, fontsize=15)

    plt.xticks(np.arange(len(ring_types)), labels = x_labels, fontweight = "bold", fontsize = 20)
    plt.yticks(np.linspace(0, 1, 6), fontsize = 20)
    plt.ylim(0.0, 1.2)
    plt.xlabel("Ring Type", fontweight = "bold", fontsize = 20, labelpad = 15)
    plt.ylabel("Density", fontsize = 20)
    plt.gca().tick_params(axis="x", pad = 120) 
    # if len(args) == 8:
    #     plt.title("Barplot Showing Proportion of Molecules with Different Ring Types with Counts in Thiols databases", fontsize = 18)
    # elif len(args) == 5:
    #     plt.title("Barplot Showing Proportion of Molecules with Different Ring Types with Counts in OE62 databases", fontsize = 18)
    # elif len(args) == 4:
    #     plt.title("Barplot Showing Proportion of Molecules with Different Ring Types with Counts in OE62+THz databases", fontsize = 18)
    # else:
    #     plt.title("Barplot Showing Proportion of Molecules with Different Ring Types with Counts in qm9 databases", fontsize = 18)

    plt.legend(loc = "upper left", mode = "expand", ncol = 5, fontsize = 16)
    plt.margins(x=0)
    plt.tight_layout()
    plt.show()
    #plt.savefig(f"/root/MChem_DGMs/analysis/Plots/DRUGS/rings.pdf", format='pdf', dpi=300)

def show_correlation(args,database, smiles_file, substructure):
    """Shows distribution of p values for a given database with and without
    the presence of a select substructure.
    """

    _, smiles_matches = find_substructure(args ,database, smiles_file, substructure)

    db = connect(f"{base_dir}/{database}")
    group_absent_ps = []
    group_present_ps = []

    for row in db.select():
        if row.SMILES in smiles_matches:
            group_present_ps.append(row.P_value)
        else:
            group_absent_ps.append(row.P_value)

    # plt.hist(group_absent_ps, alpha = 0.6, label = "Absent", density = True)
    # plt.hist(group_present_ps, alpha = 0.6, label = "Present", density = True)

    # plt.xlabel("P Value")
    # plt.ylabel("Density")
    # plt.title(f"Histogram Comparing Distribution of P Values for {substructure} Present or Absent in {database}")
    # plt.margins(x=0)
    # plt.legend()
    # plt.show()

    absent_mean = np.array(group_absent_ps).mean()
    present_mean = np.array(group_present_ps).mean()

    return absent_mean, present_mean, group_absent_ps, group_present_ps


def compare_correlations(args, substructure):
    """Plots barplots of mean p values when a substructure is present and
    absent from all databases.
    """

    absent_means = []
    present_means = []
    absent_distributions = []
    present_distributions = []
    num_datasets = len(args)

    for arg_set in args:
        absent_mean, present_mean,group_absent_ps, group_present_ps = show_correlation(args, arg_set[0], arg_set[1], substructure)
        absent_means.append(absent_mean)
        present_means.append(present_mean)
        absent_distributions.append(group_absent_ps)
        present_distributions.append(group_present_ps)

    x = 1.6*np.array(list(range(num_datasets)))
    bar_width = 0.65

    # plt.bar(x - bar_width / 2, absent_means, edgecolor = "black", align = "center", width = bar_width, label = "Absent")
    # plt.bar(x + bar_width / 2, present_means, edgecolor = "black", align = "center", width = bar_width, label = "Present")

    database_labels = [set_label(arg_set[0]) if args.index(arg_set) % 2 != 0 else set_label(arg_set[0]) for arg_set in args]
# database_labels = ["\n" + set_label(arg_set[0]) if args.index(arg_set) % 2 != 0 else set_label(arg_set[0]) for arg_set in args]

    # plt.xlabel("Database")
    # plt.ylabel("Mean P Value")
    # if len(args) == 8:
    #     plt.title(f"Barplot Showing Mean P Values of Molecules When {substructure} Present or Absent across Thiols databases")
    # elif len(args) == 5:
    #     plt.title(f"Barplot Showing Mean P Values of Molecules When {substructure} Present or Absent across OE62 databases")
    # elif len(args) == 3:
    #     plt.title(f"Barplot Showing Mean P Values of Molecules When {substructure} Present or Absent across OE62+THz databases")
    # else:
    #     plt.title(f"Barplot Showing Mean P Values of Molecules When {substructure} Present or Absent across qm9 databases")
    # plt.xticks(np.arange(num_datasets), labels = database_labels)
    # plt.margins(x = 0)
    # plt.legend()
    # plt.show()
    plt.figure(figsize=(9, 5))
    violin1=plt.violinplot(absent_distributions, x - bar_width / 2, points=60, widths=0.7,
                     showmeans=True, #showextrema=True, showmedians=True,
                     #quantiles=[0.05, 0.1, 0.8, 0.9], 
                     bw_method=0.5)#, side='low')

    violin2=plt.violinplot(present_distributions, x + bar_width / 2, points=60, widths=0.7,
                        showmeans=True, #showextrema=True, showmedians=True,
                       # quantiles=[0.05, 0.1, 0.8, 0.9], 
                        bw_method=0.5)#, side='high')
    plt.xticks(x, labels = database_labels,rotation = 40,ha='right')#,va="right")
    plt.ylabel("P")
    plt.legend(["Absent","Present"],loc=[1.05,0.4])
    plt.tight_layout()
    plt.show()



def get_functional_groups(args):
    """Creates a dictionary of functional groups present in a database.
    """


    def get_group_dict(smiles_file):
        print(smiles_file)
        with open(f"{base_dir}/{smiles_file}", "r", encoding = "cp1252") as f:
            smiles_list = json.load(f)
            group_dict = {}
            for smiles in smiles_list:
                if smiles == "INVALID":
                    continue
                else:
                    mol = Chem.MolFromSmiles(smiles)
                    fgs = identify_functional_groups(mol)
                    for group in fgs:
                        group_smiles = group.atoms
                        if group_smiles not in group_dict:
                            group_dict[group_smiles] = [smiles_list.index(smiles)]
                        else:
                            group_dict[group_smiles].append(smiles_list.index(smiles))

            return group_dict, smiles_list
            


    check = args[1]
    all_smiles_lists = []
    all_group_dicts = []

    if len(check) == 2:
        if len(args) == 2:
            arg_sets = args
        elif len(args) in (3, 8):
            arg_sets = [args[0], args[1]]
        else:
            arg_sets = [args[0], args[2]]
        for arg_set in arg_sets:
            group_dict, smiles_list = get_group_dict(arg_set[1])
            all_group_dicts.append(group_dict)
            all_smiles_lists.append(smiles_list)

        groups_list = list(all_group_dicts[0].keys())
        get_keys = lambda dictionary: list(dictionary.keys())
        other_groups = list(itertools.chain.from_iterable(list(map(get_keys, all_group_dicts[1:]))))

        for group in other_groups:
            if group not in groups_list:
                groups_list.append(group)
        
        if len(args) == 8:
            groups_file = "thiols_groups.json"
            features_file = "thiols_features.json"
        # elif len(args) == 5:
        #     groups_file = "OE62_groups.json"
        #     features_file = "OE62_features.json"
        elif len(args) == 3:
            groups_file = "OE62+THz_groups.json"
            features_file = "OE62+THz_features.json"
        else:
            groups_file = "qm9_groups.json"
            features_file = "qm9_features.json"

    else:
        smiles_file = args[1]
        print(smiles_file)
        group_dict, smiles_list = get_group_dict(smiles_file)
        groups_list = list(group_dict.keys())
        all_group_dicts.append(group_dict)
        all_smiles_lists.append(smiles_list)

        groups_file = smiles_file.replace("smiles.json", "groups.json")
        features_file = smiles_file.replace("smiles.json", "features.json")
    
    #sorted_group_dict = sorted(group_dict, key = group_dict.get, reverse = True)
    #for smiles in sorted_group_dict[:5]:
    #    print(smiles)
    #    try:
    #        mol = Chem.MolFromSmarts(smiles)
    #        image = Draw.MolToImage(mol)
    #        image.show()
    #    except:
    #        print("failed")


    def get_feature_vectors(smiles_list, group_dict):
        feature_vector_list = []
        for i, smiles in enumerate(smiles_list):
            if smiles == "INVALID":
                continue
            else:
                feature_vector = list(np.zeros(len(groups_list)))
                for j, group in enumerate(groups_list):
                    try:
                        feature_vector[j] = group_dict[group].count(i)
                    except KeyError:
                        feature_vector[j] = 0
            feature_vector_list.append(feature_vector)
        
        return feature_vector_list

    feature_vector_list_all = []
    if len(check) == 2:
        for group_dict, smiles_list in zip(all_group_dicts, all_smiles_lists):
            feature_vector_list_all.append(get_feature_vectors(smiles_list, group_dict))
    else:
        feature_vector_list_all.append(get_feature_vectors(smiles_file, group_dict))

    with open(f"/root/2192510/myenv/molecular_database_analysis/URSS-main/Functional Groups/{features_file}", "w", encoding = "cp1252") as f:
        json.dump(feature_vector_list_all, f)
    with open(f"/root/2192510/myenv/molecular_database_analysis/URSS-main/Functional Groups/{groups_file}", "w", encoding = "cp1252") as f:
        json.dump(groups_list, f)


def main():
    """The main program, where the user enters inputs to choose which 
    functions to run.
    """
    RDLogger.DisableLog('rdApp.*')

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
                args = (
                ("qm9_smiles.db", "qm9_smiles.json"),
                ("Gschnet_qm9.db", "Gschnet_qm9_smiles.json"),
                ("EDM_qm9.db","EDM_qm9_smiles.json"),
                ("GeoLDM_qm9.db", "GeoLDM_qm9_smiles.json")

                )
            else:
                args = (
               ("qm9_smiles.db", "qm9_smiles.json"),
               ("Gschnet_qm9_filtered.db", "Gschnet_qm9_filtered_smiles.json"),
               ("EDM_qm9_filtered.db","EDM_qm9_filtered_smiles.json"),
               ("GeoLDM_qm9_filtered.db","GeoLDM_qm9_filtered_smiles.json")
               )##filtered

        elif choice1 =="2":
            if filtered == 'y':##raw6
               args = (
                   ("geom_drugs.db","geom_drugs_smiles.json"),
                   ("Gschnet_drugs_filtered.db","Gschnet_drugs_filtered_smiles.json"),

                   ("EDM_drugs_filtered.db","EDM_drugs_filtered_smiles.json"),
                    ("GeoLDM_drugs_filtered.db","GeoLDM_drugs_filtered_smiles.json"),
                    ("JODO_drugs_filtered.db","JODO_drugs_filtered_smiles.json")
                )
            else:
                args = (
                    ("geom_drugs.db","geom_drugs_smiles.json"),
                    ("Gschnet_drugs_raw.db","Gschnet_drugs_raw_smiles.json"),

                    ("EDM_drugs_raw.db","EDM_drugs_raw_smiles.json"),
                    ("GeoLDM_drugs_raw.db","GeoLDM_drugs_raw_smiles.json"),
                    ("JODO_drugs_raw.db","JODO_drugs_raw_smiles.json")
                )


        elif choice1 == "3":
            if filtered == 'n':
                args = (
                ("OE62_full.db", "OE62_full_smiles.json"),
                    ("Gschnet_oe62_raw.db", "Gschnet_oe62_raw_smiles.json"),
                    ("EDM_oe62_raw.db", "EDM_oe62_raw_smiles.json")

                    )
                
            else:
                args = (
                    ("OE62_full.db", "OE62_full_smiles.json"),
                    ("Gschnet_oe62_filtered.db", "Gschnet_oe62_filtered_smiles.json"),
                    ("EDM_oe62_filtered.db", "EDM_oe62_filtered_smiles.json")
                )



        # elif choice1 == "4":
        #     args = (
        #        ("OE62_full.db", "OE62_full_ecomp.json"),
        #        ("qm9_smiles.db", "qm9_ecomp.json"),
        #         ("geom_drugs.db","geom_drugs_ecomp.json")
        #     )
            
        
        else:
            quit()

        # If the user chooses, they can write or rewrite the files containing
        # the element composition data for each database.
        # choice2 = input("Generate new smiles data (y/n)?\n")
        # if choice2 == "y":
        #     if choice1 in {"1","2","3"}:
        #         for arg_set in args:
        #             generate_smiles(arg_set[0], arg_set[1])
        #     else:
        #         generate_smiles(args[0], args[1])

        while True:
            choice3 = input(
                """Would you like to get data in terms of:
                1. How many molecules have the functional group,
                2. in terms of how many times it appears in total,
                3. miscellaneous features, e.g. no. of rings, degree of unsaturation, or
                4. how the functional group correlates with the target p value?
                5. Alternatively, identify functional groups in the database.
                Enter anything else to go back.
                """
                )

            if choice3 == "1":
                func_group = input("Enter the functional group to use (SMARTS form).\n")
                if choice1 in {"1","2","3"}:
                    for arg_set in args:
                        get_appearances(arg_set[0], arg_set[1], func_group)
                else:
                    get_appearances(args[0], args[1], func_group)

            elif choice3 == "2":
                func_group = input("Enter the functional group to use (SMARTS form).\n")
                f_group_name = input("Enter the name of this functional group.\n")

                split = False
                if choice1 in {"1","2","3"}:
                    use_sampled = input("Use sampled generated databases? (y/n)\n")
                    if use_sampled == "y":
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
                    
                        original_args = [arg_set for arg_set in args]
                        generated_samples = sample_generated(args, sample_data = data, bin_width = bin_width, split = split)
                        args = original_args
                    else:
                        generated_samples = None

                if choice1 in {"1","2","3"}:
                    create_barplot_all(args, func_group, f_group_name, generated_samples)
                else:
                    create_barplot(args[0], args[1], func_group, f_group_name)

            elif choice3 == "3":
                if choice1 in {"1","2","3"}:
                    ring_nums_all = []
                    ring_type_nums_list = []
                    aro_atoms_count_list = []
                    ring_type_counts_all = []
                    avg_ihd_list = []

                    for arg_set in args:
                        ring_nums, ring_type_nums, aro_rings_atom_counts, ring_type_counts, avg_ihd = get_misc_info(arg_set[0], arg_set[1])
                        ring_nums_all.append(ring_nums)
                        ring_type_nums_list.append(ring_type_nums)
                        aro_atoms_count_list.append(aro_rings_atom_counts)
                        ring_type_counts_all.append(ring_type_counts)
                        avg_ihd_list.append(avg_ihd)

                    plot_misc_data(args, ring_nums_all, ring_type_nums_list, aro_atoms_count_list, ring_type_counts_all, avg_ihd_list)

                else:
                    ring_nums, ring_type_nums, aro_rings_atom_counts, avg_ihd = get_misc_info(args[0], args[1])

            elif choice3 == "4":
                func_group = input("Enter the functional group to use (SMARTS form).\n")
                if choice1 in {"1","2","3"}:
                    compare_correlations(args, func_group)
                else:
                    show_correlation(args, args[0], args[1], func_group)
            
            elif choice3 == "5":
                get_functional_groups(args)
            else:
                break


if __name__ == "__main__":
    main()
