"""Provides a means of analysing the average elemental composition of
molecules or the distribution of a particular element in the ase
databases. This data can be visualised either for an individual 
or all at once.
"""
from cProfile import label
from calendar import c
import pickle 
import tkinter as tk
import json
import random
import math
from tkinter import font
import ase
from ase.db import connect
from ase.data import atomic_numbers, atomic_masses
import matplotlib as mpl
mpl.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import numpy as np
from numpy import linspace
from scipy.stats import gaussian_kde, linregress
from cycler import cycler




size=20
params = {'legend.fontsize': size,
          'figure.figsize': (7,5), #(9,5),#
          'lines.linewidth': 5,
          'axes.labelsize': size,
          'axes.titlesize': size,
          'xtick.labelsize': size,
          'ytick.labelsize': size,
          'axes.titlepad': 10,
         'axes.linewidth':1.5,
        # "axes.prop_cycle": plt.cycler("color", plt.cm.RdYlGn(np.linspace(0,1,8,endpoint=True)))
         'axes.prop_cycle': mpl.cycler(color=["k", "r", "b", "g","m"])
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

   

    elif "EDM" in database:
        label_db="EDM"

    elif "GeoLDM" in database:
        label_db="GeoLDM"
    elif "JODO" in database:
        label_db="JODO"

    if "filtered" in database:
        label_db += "_filtered"
    else: 
        label_db += ""

    
    return label_db

def generate_data(database, data_file):
    """Defining a function to generate elemental composition analysis data
    for a chosen database.
    """
    # print(sqlite3.version) 
    # print(sqlite3.sqlite_version)
    # Connecting to one of the databases.

    db = ase.db.connect(f"MChem_DGMs/analysis/Databases/{database}")

    # Initialising an empty dictionary to store the number of times each
    # element occurs in total.
    element_counts_total = {}

    #Initialising an empty list to store all data objects.
    all_element_data = []

    # Iterating through the database to get the number of atoms of each
    # element for each molecule.
    for row in db.select():
        atoms = row.toatoms()
        elements = atoms.get_chemical_symbols()
        # Initialising an empty dictionary to store the number of times each
        # element occurs in each molecule (resets every iteration)
        element_counts = {}

        # Another loop which adds the count of each element to the dictionaries.
        for i in set(elements):
            element_count = elements.count(i)
            element_counts[i] = element_count

            if i not in element_counts_total:
                element_counts_total[i] = element_count

            else:
                element_counts_total[i] += element_count

        # Adds the dictionaries containing elemental composition data for each
        # molecule to the empty list.
        all_element_data.append(element_counts)

    # Adds the dictionary containing elemental composition data for the
    # whole database to the empty list, then dumps the list to a json file.
    all_element_data.append(element_counts_total)
    f = open(f"MChem_DGMs/analysis/Databases/{data_file}", "w", encoding = "cp1252")
    json.dump(all_element_data, f)
    f.close()


def get_average_composition(data_file, use_valid = False, sample = None, split = False):
    """A function which takes the dictionaries containing the number of
    atoms of each element in each molecule, finds the percentage elemental
    composition of each and divides it by the number of molecules. The 
    result is displayed to the user.
    """

    if use_valid == False:
        f = open(f"MChem_DGMs/analysis/Databases/{data_file}", "r", encoding = "cp1252")
        ecomp_data = json.load(f)
        if sample is not None:
            summary = ecomp_data[-1]
            ecomp_data = [ecomp for i, ecomp in enumerate(ecomp_data) if i in sample]
            ecomp_data.append(summary)
        f.close()
    else:
        ecomp_data = get_valid(data_file, sample)

    if split is False:
        avg_mol_comp = {}
    else:
        avg_mol_comp = [{}, {}]
        contains_gold = []

    for molecule in ecomp_data:
        if split is True:
            if "Au" in list(molecule.keys()):
                contains_gold.append(1)
            else:
                contains_gold.append(0)
        num_atoms = sum(list(molecule.values()))
        for symbol, count in molecule.items():
            try:
                if "Au" in list(molecule.keys()) and split is True:
                    avg_mol_comp[1][symbol] += count/num_atoms
                elif split is True:
                    avg_mol_comp[0][symbol] += count/num_atoms
                else:
                    avg_mol_comp[symbol] += count/num_atoms
            except KeyError:
                if "Au" in list(molecule.keys()) and split is True:
                    avg_mol_comp[1][symbol] = count/num_atoms
                elif split is True:
                    avg_mol_comp[0][symbol] = count/num_atoms
                else:
                    avg_mol_comp[symbol] = count/num_atoms
    
    if split is False:
        for element in avg_mol_comp:
            avg_mol_comp[element] /= len(ecomp_data) - 1
    else:
        for i, avg_mol_comp_dict in enumerate(avg_mol_comp):
            for element in avg_mol_comp_dict:
                if i == 0:
                    avg_mol_comp_dict[element] /= contains_gold.count(0)
                else:
                    avg_mol_comp_dict[element] /= contains_gold.count(1)
    
    return avg_mol_comp


def get_element_counts(data_file, symbol, use_valid = False, sample = None):
    """A function which gets the count of a chosen element in each molecule."""

    if use_valid == False:
        f = open(f"MChem_DGMs/analysis/Databases/{data_file}", "r", encoding = "cp1252")
        ecomp_data = json.load(f)
        if sample is not None:
            summary = ecomp_data[-1]
            ecomp_data = [ecomp for i, ecomp in enumerate(ecomp_data) if i in sample]
            ecomp_data.append(summary)
        f.close()

    else:
        ecomp_data = get_valid(data_file, sample)

    molecule_count = len(ecomp_data) - 1
    ecount_arr = np.zeros((molecule_count,), dtype=int)

    for molecule_num in range(0, molecule_count):
        try:
            ecount_arr[molecule_num] = ecomp_data[molecule_num][symbol]
        except KeyError:
            pass

    return ecount_arr


def create_histogram(database, data_file, mode = "all"):
    """A function which plots a histogram for the distribution of a chosen
    element across molecules in a select database.
    """

    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    element = input(
        """Enter the symbol of the element you would like to create the
        histogram for (if you don't know what elements to expect, see the
        molecular composition first).
        """
        )
    
    for boolean in mode_dict[mode]:
        ecount_arr = get_element_counts(data_file, element, use_valid = boolean)
        bins = np.arange(ecount_arr.max() + 1) - 0.5
        label_db=set_label(database)
        if boolean == False:
            label = label_db# (all molecules)"
        else:
            label = label_db+"(valid only)"
        plt.hist(ecount_arr, bins = bins, edgecolor = "black", density = True, alpha = 0.6, label = label)

    plt.xlabel("Count")
    plt.ylabel("Density")
    # plt.title(f"Histogram Showing Distribution of {element} in {database}")
    plt.margins(x=0)
    plt.legend()
    plt.show()


def create_histogram_all(args, mode = "all", samples = None):
    """A function which plots a histogram showing the distribution of a
    chosen element across molecules for all databases.
    """

    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    element = input(
        """Enter the symbol of the element you would like to create the
        histogram for (if you don't know what elements to expect, see the
        molecular composition first).
        """
        )
    num_datasets = len(args)
    if mode == "compare":
        bar_width = 0.8 / (num_datasets * 2)
        shifts = np.arange(-num_datasets, num_datasets) * bar_width
    else:
        bar_width = 0.8 / num_datasets
        shifts = np.arange(-num_datasets / 2, num_datasets / 2) * bar_width + bar_width/2

    plt.figure()#figsize=(16, 12))

    
    
    for boolean in mode_dict[mode]:
        hatches = ["","","//","","//","","//"]
        for i, arg_set in enumerate(args):
            # if len(args) == 5 and i % 2 != 0:
            #     continue
            if i == 0 or samples == None:
                ecount_arr = get_element_counts(arg_set[1], element, boolean)
            else:
                ecount_arr = get_element_counts(arg_set[1], element, boolean, samples[i-1])
            bins = np.arange(ecount_arr.max() + 1) - 0.5
            shifted_bins = bins + shifts[i]
            label_db=set_label(arg_set[0])
            if boolean == False:# and mode == "compare":
                shifted_bins = bins + shifts[:num_datasets][i]
                label = label_db#f"{arg_set[0]} (all molecules)"
            else:
                label = label_db+"(valid only)"#f"{arg_set[0]} (valid only)"
                if mode == "compare":
                    shifted_bins = bins + shifts[num_datasets:][i]
                else:
                    shifted_bins = bins + shifts[i]

            data = np.array([ecount_arr, shifted_bins], dtype=object)
            np.save(f"/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/data_{element}_{label}", data)


            # n,b,patches =
            plt.hist(
                ecount_arr,
                bins = shifted_bins,
                edgecolor = "black",
                width = bar_width,
                rwidth = bar_width,
                label = label,
                density = True
                )
            # for p in patches:
            #     p.set_hatch(hatches[i])

        plt.xlabel("Count",fontsize=35)
        plt.ylabel("Density",fontsize=35)
        plt.tick_params(axis='both', labelsize=25)

    #plt.xticks(ticks=range(0, 19, 1))  # 19 for h 
    #plt.xlim(left=-0.5)

    # if len(args) == 8:
    #     plt.title(f"Histogram Showing Distribution of {element} in thiols databases")
    # elif len(args) == 2:
    #     plt.title(f"Histogram Showing Distribution of {element} in qm9 databases")
    # elif len(args) == 5:
    #     plt.title(f"Histogram Showing Distribution of {element} in OE62 databases")
    # else:
    #     plt.title(f"Histogram Showing Distribution of {element} in OE62+THz databases")
    plt.margins(x=0)
    plt.legend()
    plt.tight_layout()

    
    plt.show()

    return 


def fit_spline(database, data_file, mode = "all"):
    """A function which plots a kernel density estimate plot for the
    distribution of an element for a select database.
    """

    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    element = input(
        """Enter the symbol of the element you would like to create the 
        histogram for (if you don't know what elements to expect, see the
        molecular composition first).
        """
        )

    plt.figure()#figsize=(12, 8))

    for boolean in mode_dict[mode]:
        ecount_arr = get_element_counts(data_file, element, boolean)
        density = gaussian_kde(ecount_arr, bw_method=0.3)
        x = np.linspace(0, ecount_arr.max(), 1000)
        plt.plot(x, density(x), label=database)

    plt.xlabel("Count")
    plt.ylabel("Density")
    # plt.title(f"Kernel Density Estimate Showing Distribution of {element} in {database}")
    plt.margins(x=0)
    plt.legend()
    plt.tight_layout()
    plt.show()


def fit_spline_all(args, mode = "all", samples = None):
    """A function which plots a kernel density estimate plot for the
    distribution of an element for all databases
    """

    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    element = input("""Enter the symbol of the element you would like to create the
                    histogram for (if you don't know what elements to expect, see the
                    molecular composition first).
                    """
                    )

    plt.figure()#figsize=(12, 8))

    for boolean in mode_dict[mode]:
        for i, arg_set in enumerate(args):
            if len(args) == 5 and i % 2 != 0 and mode != "all":
                    continue
            if i == 0 or samples == None:
                ecount_arr = get_element_counts(arg_set[1], element, boolean)
            else:
                ecount_arr = get_element_counts(arg_set[1], element, boolean, samples[i-1])
            density = gaussian_kde(ecount_arr, bw_method=0.3)
            x = np.linspace(0, ecount_arr.max(), 1000)
            label_db=set_label(arg_set[0])
            if boolean == True:
                label = label_db+"(valid)"#f"{arg_set[0]} (valid only)"
            else:
                label = label_db#f"{arg_set[0]} (all molecules)"

            np.save(f"/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/kde_{element}_{label}", np.array([x, density(x)], dtype=object))
            plt.plot(x, density(x), alpha=0.7, label = label)

    plt.xlabel("Count")
    plt.ylabel("Density")
    # if len(args) == 8:
    #     plt.title(f"Kernel Density Estimate Showing Distribution of {element} in thiols databases")
    # elif len(args) == 2:
    #     plt.title(f"Kernel Density Estimate Showing Distribution of {element} in qm9 databases")
    # elif len(args) == 5:
    #     plt.title(f"Kernel Density Estimate Showing Distribution of {element} in OE62 databases")
    # else:
    #     plt.title(f"Kernel Density Estimate Showing Distribution of {element} in OE62+THz databases")
    
    plt.margins(x=0)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_hist_kde(database, data_file, mode = "all"):
    """This function plots the kde on top of a normalised histogram for a
    select database, so that it can be seen that the two distributions
    are consistent.
    """

    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    element = input(
        """Enter the symbol of the element you would like to create the
        histogram for (if you don't know what elements to expect, see the
        molecular composition first).
        """
        )

    for boolean in mode_dict[mode]:
        ecount_arr = get_element_counts(data_file, element, boolean)
        bins = np.arange(ecount_arr.max() + 1) - 0.5  # Adjust the bins to match counts
        plt.hist(ecount_arr, bins = bins, edgecolor = "black", alpha = 0.5, density = True)

        density = gaussian_kde(ecount_arr, bw_method=0.3)
        x = np.linspace(0, ecount_arr.max(), 1000)
        plt.plot(x, density(x), alpha=0.7)

    plt.xlabel("Count")
    plt.ylabel("Density")
    # plt.title(f"Histogram and KDE Showing Distribution of {element} in {database}")
    plt.margins(x=0)
    plt.tight_layout()
    plt.show()


def plot_kde_hist_all(args, mode = "all", samples = None):
    """This function plots the kde on top of a normalised histogram for all
    databases, so that it can be seen that the two distributions are
    consistent.
    """

    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    element = input(
        """Enter the symbol of the element you would like
        to create the histogram for (if you don't know what elements to
        expect, see the molecular composition first).
        """
        )
    num_datasets = len(args)
    if mode == "compare":
        bar_width = 0.8 / (num_datasets * 2)
        shifts = np.arange(-num_datasets, num_datasets) * bar_width
    else:
        bar_width = 0.8 / num_datasets
        shifts = np.arange(-num_datasets / 2, num_datasets / 2) * bar_width

    plt.figure()#figsize=(16, 12))

    for boolean in mode_dict[mode]:
        for i, arg_set in enumerate(args):
            if len(args) == 5 and i % 2 != 0 and mode != "all":
                continue
            # Plotting histogram
            if i == 0 or samples == None:
                ecount_arr = get_element_counts(arg_set[1], element, boolean)
            else:
                ecount_arr = get_element_counts(arg_set[1], element, boolean, samples[i-1])
            bins = np.arange(ecount_arr.max() + 1) - 0.5
            label_db=set_label(arg_set[0])
            if boolean == False:
                shifted_bins = bins + shifts[:int(num_datasets)][i]
                label = label_db#f"{arg_set[0]} (all molecules)"
            else:
                label = label_db+"(valid)"# f"{arg_set[0]} (valid only)"
                if mode == "compare":
                    shifted_bins = bins + shifts[int(num_datasets):][i]
                else:
                    shifted_bins = bins + shifts[i]
            plt.hist(
                ecount_arr,
                bins = shifted_bins,
                edgecolor = "black",
                width = bar_width,
                rwidth = bar_width,
                label = label,
                density = True
                )

            # Plotting KDE
            density = gaussian_kde(ecount_arr, bw_method=0.3)
            x = np.linspace(0, ecount_arr.max(), 1000)
            plt.plot(x, density(x), label=arg_set[0], alpha=0.7)

    plt.xlabel("Count")
    plt.ylabel("Density")
    # if len(args) == 8:
    #     plt.title(f"Histogram and KDE Showing Distribution of {element} in thiols databases")
    # elif len(args) == 2:
    #     plt.title(f"Histogram and KDE Showing Distribution of {element} in qm9 databases")
    # elif len(args) == 5:
    #     plt.title(f"Histogram and KDE Showing Distribution of {element} in OE62 databases")
    # else:
    #     plt.title(f"Histogram and KDE Showing Distribution of {element} in OE62+THz databases")

    plt.margins(x=0)
    plt.legend()
    plt.tight_layout()
    plt.show()


def get_average_composition_num(data_file, use_valid, sample = None, split = False):
    """A function which converts the symbol keys in the average molecular
    composition dictionary to atomic numbers.
    """

    avg_mol_comp = get_average_composition(data_file, use_valid, sample, split)
    if split is True:
        avg_mol_comp_num = [{}, {}]
        for i, avg_mol_comp_dict in enumerate(avg_mol_comp):
            for element, count in avg_mol_comp_dict.items():
                atomic_num = atomic_numbers[element]
                avg_mol_comp_num[i][atomic_num] = count
    else:
        avg_mol_comp_num = {}
        for element, count in avg_mol_comp.items():
            atomic_num = atomic_numbers[element]
            avg_mol_comp_num[atomic_num] = count
    
    return avg_mol_comp_num


def create_bar_comp(database, data_file, mode = "all"):
    """A function which generates a barplot showing average atom counts per
    molecule in a select database. This is done in terms of atomic number
    so it is easier to compare average size of molecules across databases.
    """

    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    
    for boolean in mode_dict[mode]:
        avg_mol_comp_num = get_average_composition_num(data_file, boolean)
        atomic_nums = list(map(str, sorted(list(avg_mol_comp_num.keys()))))
        average_count = list(dict(sorted(avg_mol_comp_num.items())).values())
        plt.bar(atomic_nums, average_count, edgecolor='black', alpha = 0.6,width=0.03)

    plt.xlabel('Atomic Number')
    plt.ylabel('Average Proportion per Molecule')
    # plt.title(f'Barplot Showing Average Percent Composition Per Molecule in {database}')
    plt.margins(x=0)
    plt.tight_layout()
    plt.show()


def create_bar_all(args, mode = "all", samples = None, split = False):
    """A function which generates a barplot showing average atom percentage
    per molecule in all databases. This is done in terms of atomic number
    so it is easier to compare average size of molecules between
    databases.
    """

    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    if len(args) == 5 and mode != "all":
        num_datasets = 3
    else:
        if split == True and len(args) == 3:
            args[0][0] = "OE62"
            args.insert(1, ("THz", args[0][1]))
            args[2][0] = "Generated\n(transform)"
            args[3][0] = "Generated\n(no transform)"
        num_datasets = len(args)
    all_elements = []
    avg_mol_comp_list = []
    atomic_nums_list = []

    # The purpose of this loop is primarily to get a list of all of the
    # elements across all of the databases. Lists of other data like
    # the dictionary of elements to average molecular composition and lists
    # of atomic numbers present in each database are obtained along the way.
    for boolean in mode_dict[mode]:
        for i, arg_set in enumerate(args):
            # if (len(args) == 5 and i % 2 != 0 and mode != "all") or (len(args) == 4 and i == 1):
            #     continue
            if i == 0 or samples == None:
                avg_mol_comp_num = get_average_composition_num(arg_set[1], boolean, split = split)
            else:
                # if len(args) == 4:
                #     i -= 1
                avg_mol_comp_num = get_average_composition_num(arg_set[1], boolean, samples[i-1], split = split)

            if split is True:
                for avg_mol_comp_num_dict in avg_mol_comp_num:
                    atomic_nums_list.append(list(map(str, sorted(list(avg_mol_comp_num_dict.keys())))))
                    avg_mol_comp_list.append(avg_mol_comp_num_dict)
                    if set(atomic_nums_list[-1]) != set(all_elements):
                        all_elements += list(set(atomic_nums_list[-1]).difference(set(all_elements)))
                split = False
            else:
                atomic_nums_list.append(list(map(str, sorted(list(avg_mol_comp_num.keys())))))
                avg_mol_comp_list.append(avg_mol_comp_num)
                if set(atomic_nums_list[-1]) != set(all_elements):
                    all_elements += list(set(atomic_nums_list[-1]).difference(set(all_elements)))

        # The list of all atomic numbers is sorted as though they were integer values.
        all_elements = sorted(all_elements, key = int)

        # If there are any elements which exist in one of the databases but
        # not another, a key is created and set to 0 where it is absent.
        an_index = 0
        for atomic_nums in atomic_nums_list:
            for element in all_elements:
                if element not in atomic_nums:
                    avg_mol_comp_list[an_index][int(element)] = 0
            an_index += 1

        # x and y values obtained
        average_counts = list(map(lambda x:list(dict(sorted(x.items())).values()), avg_mol_comp_list))
        num_av_counts = np.array(average_counts)
        x = np.arange(len(all_elements))   
        # np.save("my_array",num_av_counts)
        # np.save("x_arr",x)


        # Adjusting the positions of the different plots
        # so they are all visible and changing their size so they all fit.
        plt.figure(figsize=(9,5.5))
        # N=9
        # plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.RdYlGn(np.linspace(0,1,N)))

        bar_width = 0.1#0.3
        if mode == "compare":
            shifts = np.arange(-num_datasets, num_datasets) * bar_width + bar_width/2
        else:
            shifts = np.arange(-num_datasets / 2, num_datasets / 2) * bar_width + bar_width/2
        if not(mode == "compare" and boolean == True):
            shifted_xs = [x + shifts[:num_datasets][i] for i in range(num_datasets)]
        else:
            shifted_xs = [x + shifts[num_datasets:][i] for i in range(num_datasets)]

        # Finally, the plot is made.
        for i in range(num_datasets):
            j = i
            label_db=set_label(args[i][0])
            if boolean == True and mode == "compare":
                j += num_datasets
                label = label_db #f"{args[i][0]} (valid only)"
            else:
                label = label_db #f"{args[i][0]} (all molecules)"
            plt.bar(
                shifted_xs[i],
                average_counts[j],
                bar_width,
                label = label,
                edgecolor = "black"
                )

    plt.xlabel('Atomic Number',fontsize=35)
    plt.ylabel('Average Proportion per Molecule',fontsize=30)
    # if len(args) == 8:
    #     plt.title('Barplot Showing Average Percent Composition Per Molecule in Thiols Databases')
    # elif len(args) == 2:
    #     plt.title('Barplot Showing Average Percent Composition Per Molecule in qm9 Databases')
    # elif len(args) == 5:
    #     plt.title('Barplot Showing Average Percent Composition Per Molecule in OE62 Databases')
    # else:
    #     plt.title('Barplot Showing Average Percent Composition Per Molecule in OE62+THz Databases')

    plt.xticks(np.arange(len(all_elements)), labels = all_elements, fontsize=35)
    plt.yticks(fontsize=35)
    plt.margins(x=0)
    plt.legend(fontsize=35)#(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.show()
    return x, num_av_counts


def get_atoms_data(data_file, use_valid, sample = None, split = False, homo_lumo=False):
    """A function which gets the number of atoms in, and molecular weight 
    of, each molecule in a database.
    """

    if use_valid == False:
        f = open(f"MChem_DGMs/analysis/Databases/{data_file}", "r", encoding = "cp1252")
        ecomp_data = json.load(f)
        if sample is not None:
            summary = ecomp_data[-1]
            ecomp_data = [ecomp for i, ecomp in enumerate(ecomp_data) if i in sample]
            ecomp_data.append(summary)
        f.close()
        # loading in homo-lumo gap predictions for oe62
        if homo_lumo:
            if "of_" in data_file:
                hl=np.load("/home/zkoczor//Work2/cgschnet-thiols/gschnet/OE62/homo_lumo/original_training/homo_lumo.npz",allow_pickle=True)["HL"][()]
    
            if "o2_" in data_file:
                hl=np.load("/home/zkoczor//Work2/cgschnet-thiols/gschnet/OE62/homo_lumo/OE62_2_filtered_generated/homo_lumo.npz",allow_pickle=True)["HL"][()]
    # print(min(hl),max(hl)) 
    # hl_train=np.load("/home/zkoczor//Work2/cgschnet-thiols/gschnet/OE62/homo_lumo/original_training/homo_lumo.npz",allow_pickle=True)["HL"][()]

            if sample is not None:
                hl = [ecomp for i, ecomp in enumerate(hl) if i in sample]
    else:
        ecomp_data = get_valid(data_file, sample)

    if split == True:
        no_of_atoms = [[], []]
        weights = [[], []]
        for molecule in ecomp_data[:-1]:
            if "Au" in molecule:
                no_of_atoms[1].append(sum(list(molecule.values())))
                molecular_weight = 0
                for symbol, ecount in molecule.items():
                    molecular_weight += atomic_masses[atomic_numbers[symbol]] * ecount
                weights[1].append(molecular_weight)
            else:
                no_of_atoms[0].append(sum(list(molecule.values())))
                molecular_weight = 0
                for symbol, ecount in molecule.items():
                    molecular_weight += atomic_masses[atomic_numbers[symbol]] * ecount
                weights[0].append(molecular_weight)

    else:    
        no_of_atoms = [sum(list(molecule.values())) for molecule in ecomp_data[:-1]]
        weights = []
        for molecule in ecomp_data[:-1]:
            molecular_weight = 0
            for symbol, ecount in molecule.items():
                molecular_weight += atomic_masses[atomic_numbers[symbol]] * ecount
            weights.append(molecular_weight)

    if homo_lumo:
        return no_of_atoms, weights, hl
    return no_of_atoms, weights


def ellipse(semi_major_ax, semi_minor_ax, rotation_angle, x_origin, y_origin, points = 100):
    """This is a function which is used to draw the elipses which will 
    indicate where 1, 2 and 3 * standard deviations are on the density plot.
    """

    cos_ra, sin_ra = np.cos(rotation_angle), np.sin(rotation_angle)
    theta = linspace(0, 2 * np.pi, points)
    X = semi_major_ax * np.cos(theta) * cos_ra - sin_ra * semi_minor_ax * np.sin(theta) + x_origin
    Y = semi_major_ax * np.cos(theta) * sin_ra + cos_ra * semi_minor_ax * np.sin(theta) + y_origin
    return X, Y


def plot_2D_density_atom_data(database, x, y, x_lims, y_lims):
    """Creates a 2D density plot of number of atoms vs molecular weight to
    better visualise the relationship. Additionally, the histograms are 
    shown along the x and y axes to view the individual distributions.
    """

    xlabel = "Number of Atoms"
    ylabel = 'Molecular Weight'

    # Set dimentsions of plots
    left, width = 0.17, 0.6
    bottom, height = 0.12, 0.55
    bottom_h = left + width - 0.08
    left_h = left + width + 0.02

    rect_temperature = [left, bottom, width, height]
    rect_hist_x = [left + 0.033, bottom_h, width - 0.065, 0.15]
    rect_hist_y = [left_h, bottom, 0.15, height]

    fig = plt.figure(1)#, figsize=(9.5,9))

    # Create 3 separate plots
    ax_temperature = plt.axes(rect_temperature)
    ax_hist_x = plt.axes(rect_hist_x)
    ax_hist_y = plt.axes(rect_hist_y)

    nullfmt = NullFormatter()
    ax_hist_x.xaxis.set_major_formatter(nullfmt)
    ax_hist_y.yaxis.set_major_formatter(nullfmt)

    x_min = min(x_lims)
    x_max = max(x_lims)
    y_min = min(y_lims)
    y_max = max(y_lims)

    num_x_bins = num_y_bins = x_max - x_min + 1

    # Create evenly spaced bins (there will be as many bins as there are 
    # numbers between maximum and minimum x inclusive) 
    x_bins = linspace(start = x_min, stop = x_max, num = num_x_bins)
    y_bins = linspace(start = y_min, stop = y_max, num = num_y_bins)
    aspect_ratio = 1.0*(x_max - 0)/(1 * y_max - 0)

    num_samples, x_edges, y_edges = np.histogram2d(y, x, bins = (y_bins, x_bins))
        
    # Plot the density data
    cax = ax_temperature.imshow(num_samples, extent=[x_min, x_max, y_min, 
                                                          y_max],
        interpolation = "nearest", origin = "lower", aspect = aspect_ratio)
    
    cax.set_cmap("BuPu")

    contour_color = 'black'
    x_center = np.mean(x)
    y_center = np.mean(y)
    semi_major_ax = np.std(x)
    semi_minor_ax = np.std(y)
    rotation_angle = 0

    #Plot ellipses to show 1, 2 and 3 standard deviations from the mean
    X, Y = ellipse(semi_major_ax, semi_minor_ax, rotation_angle, x_center, y_center)
    ax_temperature.plot(X, Y, "k:", color = contour_color, ms = 1, linewidth = 2)
    ax_temperature.annotate("$1\\sigma$", xy = (X[15], Y[15]), xycoords = "data", xytext = (10, 10),
                            textcoords = "offset points", horizontalalignment = "right",
                            verticalalignment = "bottom", color = contour_color, fontsize = 18
                            )
        
    X, Y = ellipse(2 * semi_major_ax, 2 * semi_minor_ax, rotation_angle, x_center, y_center)
    ax_temperature.plot(X, Y, "k:", color = contour_color, ms = 1, linewidth=2)
    ax_temperature.annotate("$2\\sigma$", xy = (X[15], Y[15]), xycoords = "data", xytext = (10, 10),
                            textcoords = "offset points", horizontalalignment = "right",
                            verticalalignment = "bottom", color = contour_color, fontsize = 18
                            )
        
    X, Y = ellipse(3 * semi_major_ax, 3 * semi_minor_ax, rotation_angle, x_center, y_center)
    ax_temperature.plot(X,Y,"k:", color = contour_color, ms=1,linewidth=2.0)
    ax_temperature.annotate("$3\\sigma$", xy = (X[15], Y[15]), xycoords = "data", xytext = (10, 10),
                            textcoords = "offset points", horizontalalignment = "right",
                            verticalalignment = "bottom", color = contour_color, fontsize = 18
                            )

    slope, intercept, r_value, p_value, std_err = linregress(x, y)
        
    ax_temperature.set_xlabel(xlabel, fontsize = 18)
    ax_temperature.set_ylabel(ylabel, fontsize = 18)

    ax_temperature.set_xlim(x_lims)
    ax_temperature.set_ylim(y_lims)
        
    #Plot histograms
    x_bins = np.arange(x_min, x_max + 2)
    y_bins = np.arange(y_min, y_max + 2)
        
    ax_hist_x.hist(x, bins = x_bins, color = "grey", density = True)
    ax_hist_y.hist(y, bins = y_bins, orientation = "horizontal", color = "grey", density = True)
    ax_hist_x.axvline(x_center, color = "k", linestyle = "dashed", linewidth = 1)
    ax_hist_y.axhline(y_center, color = "k", linestyle = "dashed", linewidth = 1)
        
    ax_hist_x.set_xlim(x_lims)
    ax_hist_y.set_ylim(y_lims)

    tick_labels = ax_temperature.get_xticklabels()
    for label in tick_labels:
        label.set_fontsize(18)
        
    tick_labels = ax_temperature.get_yticklabels()
    for label in tick_labels:
        label.set_fontsize(18)
    
    tick_labels = ax_hist_x.get_yticklabels()
    for label in tick_labels:
        label.set_fontsize(18)
    
    tick_labels = ax_hist_y.get_xticklabels()
    for label in tick_labels:
        label.set_fontsize(18)
    
    for tick in ax_hist_y.xaxis.get_major_ticks()[0::2]:
        tick.set_pad(20)
    
    #Show the plot
    plt.draw()
    # plt.title(f"2D Density Plot Showing the Relationship Between Number of\nAtoms and Molecular Weight in {database}", y = 1.36, x = -1.8, fontdict = {"fontsize": 18})

    print(f"slope: {slope}\nintercept: {intercept}\nr value: {r_value}\nstandard error: {std_err}")
    plt.tight_layout()
    plt.show()


def plot_atoms_data(database, data_file, mode = "all"):
    """Plots distributions for number of atoms and molecular weights in a
    database.
    """

    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    for boolean in mode_dict[mode]:
        no_of_atoms, weights = get_atoms_data(data_file, boolean)
        plt.hist(no_of_atoms, edgecolor = "black", density = True)

    plt.xlabel("Count")
    plt.ylabel("Density")
    # plt.title(f"Histogram Showing Distribution of Number of Atoms in {database}")
    plt.margins(x=0)
    plt.tight_layout()
    plt.show()

    for boolean in mode_dict[mode]:
        no_of_atoms, weights = get_atoms_data(data_file, boolean)
        plt.hist(weights, edgecolor = "black", density = True)
    
    plt.xlabel("Count")
    plt.ylabel("Density")
    # plt.title(f"Histogram Showing Distribution of Molecular Weights in {database}")
    plt.margins(x=0)
    plt.tight_layout()
    plt.show()

    for boolean in mode_dict[mode]:
        no_of_atoms, weights = get_atoms_data(data_file, boolean)

        x_lims = [max(no_of_atoms), min(no_of_atoms)]
        y_lims = [max(weights), min(weights)]

        plot_2D_density_atom_data(database, no_of_atoms, weights, x_lims, y_lims)

def plot_atoms_data_all(args, mode = "all", samples = None, split = False):
    """Plots distributions for number of atoms and molecular weights in all
    databases.
    """
    homo_lumo=False
    split_dbs = False
    if split == True and len(args) == 3:
        args[0][0] = "OE62"
        args.insert(1, ("THz", args[0][1]))
        args[2][0] = "Generated\n(transform)"
        args[3][0] = "Generated\n(no transform)"
        split_dbs = True
    num_datasets = len(args)
    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    plt.figure(figsize=(18, 8))
    no_of_atoms_list = []
    weights_list = []
    hl_list = []
    x_lims = [np.inf, 0]
    y_lims = [np.inf, 0]

    for boolean in mode_dict[mode]:
        for i, arg_set in enumerate(args):
            if (len(args) == 5 and i % 2 != 0 and mode != "all" and split_dbs == False) or (split_dbs == True and i == 1):
                split = False
                continue
            if ((i == 0 or samples == None) and split == False):
                if homo_lumo:
                    no_of_atoms, weights,hl = get_atoms_data(arg_set[1], boolean, split = split,homo_lumo=homo_lumo)
                    hl_list.append(hl)                
                else:
                    no_of_atoms, weights = get_atoms_data(arg_set[1], boolean, split = split)
                
                no_of_atoms_list.append(no_of_atoms)
                weights_list.append(weights)
                
            elif split == True and i == 0:
                no_of_atoms, weights = get_atoms_data(arg_set[1], boolean, split = split)
                for i in range(0, 2):
                    if samples != None:
                        no_of_atoms[i] = [atoms for j, atoms in enumerate(no_of_atoms[i]) if j in samples[i]]
                        weights[i] = [weight for j, weight in enumerate(weights[i]) if j in samples[i]]
                    no_of_atoms_list.append(no_of_atoms[i])
                    weights_list.append(weights[i])

                    if max(no_of_atoms[i]) > x_lims[1]:
                        x_lims[1] = max(no_of_atoms[i])
                    if max(weights[i]) > y_lims[1]:
                        y_lims[1] = max(weights[i])
                    if min(no_of_atoms[i]) < x_lims[0]:
                        x_lims[0] = min(no_of_atoms[i])
                    if min(weights[i]) < y_lims[0]:
                        y_lims[0] = min(weights[i])
            else:
                if split_dbs == True:
                    no_of_atoms, weights = get_atoms_data(arg_set[1], boolean, sample = samples[i], split = split)
                elif homo_lumo:
                    no_of_atoms, weights, hl = get_atoms_data(arg_set[1], boolean, sample = samples[i - 1], split = split,homo_lumo=homo_lumo)
                    hl_list.append(hl)
                    print("HOMO LUMO ", np.shape(weights_list), np.shape(hl_list))
                else:
                    no_of_atoms, weights = get_atoms_data(arg_set[1], boolean, sample = samples[i - 1], split = split)

                no_of_atoms_list.append(no_of_atoms)
                weights_list.append(weights)
                
                
                
                if max(no_of_atoms) > x_lims[1]:
                    x_lims[1] = max(no_of_atoms)
                if max(weights) > y_lims[1]:
                    y_lims[1] = max(weights)
                if min(no_of_atoms) < x_lims[0]:
                    x_lims[0] = min(no_of_atoms)
                if min(weights) < y_lims[0]:
                    y_lims[0] = min(weights)

        data = np.array([no_of_atoms_list],dtype=object)
        np.save(f"/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/no_of_atoms_unfiltered", data)

        for i, arg_set in enumerate(args):
            if (len(args) == 5 and i % 2 != 0 and mode != "all"):
                continue
            density = gaussian_kde(no_of_atoms_list[i], bw_method=0.3)
            x = np.linspace(0, max(no_of_atoms_list[i]), 1000)
            label_db=set_label(arg_set[0])
            if boolean == True:
                label = label_db#f"{arg_set[0]} (valid only)"
            else:
                label = label_db #f"{arg_set[0]} (all molecules)"
            plt.plot(x, density(x), label=label, alpha=0.7)
            print("mean value: ", x[np.argmax(density(x))], np.mean(no_of_atoms_list[i]))

    plt.xlabel("Count", fontsize = 40)
    plt.ylabel("Density", fontsize = 40)
    # if len(args) == 8:
    #     plt.title("Kernel Density Estimate Showing Distribution of Number of Atoms in thiols databases")
    # elif len(args) == 2:
    #     plt.title("Kernel Density Estimate Showing Distribution of Number of Atoms in qm9 databases")
    # elif len(args) == 5:
    #     plt.title("Kernel Density Estimate Showing Distribution of Number of Atoms in OE62 databases")
    # else:
    #     plt.title("Kernel Density Estimate Showing Distribution of Number of Atoms in OE62+THz databases")
    
    plt.margins(x=0)
    plt.xlim(0,70)#loc=[1.05,0])
    plt.legend(fontsize=30)#(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.tight_layout()
    plt.show()

    if len(args) == 5 and mode != "all":
        args = [args[i] for i in range(0, len(args)) if i % 2 == 0]
    if mode == "compare":
        args = args * 2

    plt.figure(figsize=(18, 8))

    for i, arg_set in enumerate(args):
        density = gaussian_kde(weights_list[i], bw_method=0.3)
        x = np.linspace(0, max(weights_list[i]), 1000)
        label_db=set_label(arg_set[0])
        if boolean == True:
            label = label_db #f"{arg_set[0]} (valid only)"
        else:
            label = label_db #f"{arg_set[0]} (all molecules)"
        plt.plot(x, density(x), label=label, alpha=0.7)

    plt.xlabel("Weight (g/mol)", fontsize = 40)
    plt.ylabel("Density", fontsize = 40)
    # if num_datasets == 8:
    #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in thiols databases")
    # elif num_datasets == 2:
    #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in qm9 databases")
    # elif num_datasets == 5:
    #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in OE62 databases")
    # else:
        # plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in OE62+THz databases")
    plt.xlim(70,1000)#170)#
    plt.margins(x=0)
    plt.legend(fontsize=30)#(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.tight_layout()
    plt.savefig(f"/root/MChem_DGMs/analysis/Plots/DRUGS/weights.pdf",format='pdf',dpi=300)
    plt.show()

   # plot HOMO-LUMO gap prediction
    if homo_lumo:
        for i, arg_set in enumerate(args):
            density = gaussian_kde(hl_list[i], bw_method=0.3)
            x = np.linspace(0, max(hl_list[i]), 1000)
            label_db=set_label(arg_set[0])
            if boolean == True:
                label = label_db #f"{arg_set[0]} (valid only)"
            else:
                label = label_db #f"{arg_set[0]} (all molecules)"
            plt.plot(x, density(x), label=label, alpha=0.7)

        plt.xlabel("HOMO-LUMO gap /eV")
        plt.ylabel("Density")
        # if num_datasets == 8:
        #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in thiols databases")
        # elif num_datasets == 2:
        #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in qm9 databases")
        # elif num_datasets == 5:
        #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in OE62 databases")
        # else:
            # plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in OE62+THz databases")
        plt.xlim(0,17.5)#170)#
        plt.margins(x=0)
        plt.legend()
        plt.tight_layout()
        plt.show()


    # Creating a 2D density plot with accompanying histograms on respective
    # axes for each database.
    for i, arg_set in enumerate(args):
        arg_set = list(arg_set)
        label_db=set_label(arg_set[0])
        if (mode == "compare" and i >= len(args)/2) or mode == "valid":
            arg_set[0] = label_db #f"{arg_set[0]} (valid only)"
        else:
            arg_set[0] = label_db #f"{arg_set[0]} (all molecules)"
        # plot_2D_density_atom_data(arg_set[0], no_of_atoms_list[i], weights_list[i], x_lims, y_lims)

def plot_atoms_data_all2(args, mode = "all", samples = None, split = False):
    """Plots distributions for number of atoms and molecular weights in all
    databases.
    """

    split_dbs = False
    if split == True and len(args) == 3:
        args[0][0] = "OE62"
        args.insert(1, ("THz", args[0][1]))
        args[2][0] = "Generated\n(transform)"
        args[3][0] = "Generated\n(no transform)"
        split_dbs = True
    num_datasets = len(args)
    mode_dict = {"compare": [False, True], "valid": [True], "all": [False]}

    plt.figure()#figsize=(12, 8))
    no_of_atoms_list = []
    weights_list = []
    x_lims = [np.inf, 0]
    y_lims = [np.inf, 0]

    for boolean in mode_dict[mode]:
        for i, arg_set in enumerate(args):
            if (len(args) == 5 and i % 2 != 0 and mode != "all" and split_dbs == False) or (split_dbs == True and i == 1):
                split = False
                continue
            if ((i == 0 or samples == None) and split == False):
                no_of_atoms, weights = get_atoms_data(arg_set[1], boolean, split = split)
            elif split == True and samples != None and i == 0:
                no_of_atoms, weights = get_atoms_data(arg_set[1], boolean, split = split)
                for i in range(0, 2):
                    no_of_atoms[i] = [atoms for j, atoms in enumerate(no_of_atoms[i]) if j in samples[i]]
                    weights[i] = [weight for j, weight in enumerate(weights[i]) if j in samples[i]]
                    no_of_atoms_list.append(no_of_atoms[i])
                    weights_list.append(weights[i])

                    if max(no_of_atoms[i]) > x_lims[1]:
                        x_lims[1] = max(no_of_atoms[i])
                    if max(weights[i]) > y_lims[1]:
                        y_lims[1] = max(weights[i])
                    if min(no_of_atoms[i]) < x_lims[0]:
                        x_lims[0] = min(no_of_atoms[i])
                    if min(weights[i]) < y_lims[0]:
                        y_lims[0] = min(weights[i])
            else:
                if split_dbs == True:
                    no_of_atoms, weights = get_atoms_data(arg_set[1], boolean, sample = samples[i], split = split)
                else:
                    no_of_atoms, weights = get_atoms_data(arg_set[1], boolean, sample = samples[i], split = split)
                no_of_atoms_list.append(no_of_atoms)
                weights_list.append(weights)

                if max(no_of_atoms) > x_lims[1]:
                    x_lims[1] = max(no_of_atoms)
                if max(weights) > y_lims[1]:
                    y_lims[1] = max(weights)
                if min(no_of_atoms) < x_lims[0]:
                    x_lims[0] = min(no_of_atoms)
                if min(weights) < y_lims[0]:
                    y_lims[0] = min(weights)

        for i, arg_set in enumerate(args):
            if (len(args) == 5 and i % 2 != 0 and mode != "all"):
                continue
            density = gaussian_kde(no_of_atoms_list[i], bw_method=0.3)
            x = np.linspace(0, max(no_of_atoms_list[i]), 1000)
            label_db=set_label(arg_set[0])
            if boolean == True:
                label = label_db#f"{arg_set[0]} (valid only)"
            else:
                label = label_db#f"{arg_set[0]} (all molecules)"
            plt.plot(x, density(x), label=label, alpha=0.7)

    plt.xlabel("Count")
    plt.ylabel("Density")
    # if len(args) == 8:
    #     plt.title("Kernel Density Estimate Showing Distribution of Number of Atoms in thiols databases")
    # elif len(args) == 2:
    #     plt.title("Kernel Density Estimate Showing Distribution of Number of Atoms in qm9 databases")
    # elif len(args) == 5:
    #     plt.title("Kernel Density Estimate Showing Distribution of Number of Atoms in OE62 databases")
    # else:
    #     plt.title("Kernel Density Estimate Showing Distribution of Number of Atoms in OE62+THz databases")
    
    plt.margins(x=0)
    plt.legend()
    plt.tight_layout()
    plt.show()

    if len(args) == 5 and mode != "all":
        args = [args[i] for i in range(0, len(args)) if i % 2 == 0]
    if mode == "compare":
        args = args * 2

    for i, arg_set in enumerate(args):
        density = gaussian_kde(weights_list[i], bw_method=0.3)
        x = np.linspace(0, max(weights_list[i]), 1000)
        label_db=set_label(arg_set[0])
        if boolean == True:
            label = label_db #f"{arg_set[0]} (valid only)"
        else:
            label = label_db #f"{arg_set[0]} (all molecules)"
        plt.plot(x, density(x), label=label, alpha=0.7)

    plt.xlabel("Weight")
    plt.ylabel("Density")
    # if num_datasets == 8:
    #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in thiols databases")
    # elif num_datasets == 2:
    #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in qm9 databases")
    # elif num_datasets == 5:
    #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in OE62 databases")
    # else:
    #     plt.title("Kernel Density Estimate Showing Distribution of Molecular Weights in OE62+THz databases")
    plt.margins(x=0)
    plt.legend()
    plt.tight_layout()
    plt.show()

    # Creating a 2D density plot with accompanying histograms on respective
    # axes for each database.
    for i, arg_set in enumerate(args):
        arg_set = list(arg_set)
        label_db=set_label(arg_set[0])
        if (mode == "compare" and i >= len(args)/2) or mode == "valid":
            arg_set[0] = label_db#f"{arg_set[0]} (valid only)"
        else:
            arg_set[0] = label_db#f"{arg_set[0]} (all molecules)"
        plot_2D_density_atom_data(arg_set[0], no_of_atoms_list[i], weights_list[i], x_lims, y_lims)
 

def get_valid(data_file, sample = None):
    """Gets a list of elemental compositions for valid molecules"""

    f = open(f"MChem_DGMs/analysis/Databases/{data_file}", "r", encoding = "cp1252")
    ecomp_data = json.load(f)
    # print(ecomp_data)
    if sample is not None:
        summary = ecomp_data[-1]
        ecomp_data = [ecomp for i, ecomp in enumerate(ecomp_data) if i in sample]
        ecomp_data.append(summary)

    smiles_file = data_file.replace("ecomp.json", "smiles.json")
    f = open(f"MChem_DGMs/analysis/Databases/{smiles_file}", "r", encoding = "cp1252")
    smiles_data = json.load(f)
    # print("SMILES ",smiles_data)
    if sample is not None:
        smiles_data = [smiles for i, smiles in enumerate(smiles_data) if i in sample]

    valid = []
    for i, j in enumerate(ecomp_data[:-1]):
        if smiles_data[i] != "INVALID":
            valid.append(j)
    
    valid.append(ecomp_data[-1])
    
    return valid


def sample_generated(args, plot = False, sample_data = "Molecular Weights", bin_width = 10, min_sample_size = 100, split = False, use_valid = False):
    """A function to generate subsets of generated databases by
    selecting a sample such that molecular weight distributions match the
    original database.
    """

    generated_samples = []
    # if split == True and len(args) == 3:
    #     args.insert(1, ["OE62", args[0][1]])
    #     args.insert(2, ["THz", args[0][1]])
    #     args[3][0] = "Generated\n(transform)"
    #     args[4][0] = "Generated\n(no transform)"

    for arg_set in args:
        print(arg_set[1])
        if args.index(arg_set) == 2 and split == True:
            split = False
            continue
        if sample_data == "Molecular Weights":
            if args.index(arg_set) == 0 and split == True:
                
                _, data = get_atoms_data(arg_set[1].replace("smiles", "ecomp"), use_valid = use_valid)
            else:
                _, data = get_atoms_data(arg_set[1].replace("smiles", "ecomp"), use_valid = use_valid, split = split)
        else:
            if args.index(arg_set) == 0 and split == True:
                data, _ = get_atoms_data(arg_set[1].replace("smiles", "ecomp"), use_valid = use_valid)
            else:
                data, _ = get_atoms_data(arg_set[1].replace("smiles", "ecomp"), use_valid = use_valid, split = split)
        if len(data) != 2:
            data_list = [data]
        else:
            data_list = data
        # Runs only for the training database - this gets the densities that the
        # generated databases should try to match.
        for data in data_list:
            if args.index(arg_set) == 0:
                num_bins = math.ceil(max(data) / bin_width)
                bins = bin_width * np.arange(num_bins + 1) - 0.5
                heights, _, _ = plt.hist(data, bins = bins, edgecolor = "black", density = True)
                train_densities = heights * bin_width
                plt.close()
            # Runs only for generated databases.
            else:
                #Begin by creating a dictionary to group weights into bins.
                bin_indices = np.digitize(data, bins) - 1
                grouped_data = {i: [] for i in range(0, len(train_densities))}
                data_indices = [i for i in range(0, len(data))]
                for idx, bin_idx in zip(data_indices, bin_indices):
                    try:
                        grouped_data[bin_idx].append(idx)
                    except KeyError:
                        continue
                
                # The bin with the biggest percentage decrease between the generated and
                # training databases is found.
                gen_densities = {i: len(group)/len(data) for i, group in grouped_data.items()}
                generated_sample = []
                previous = []
                while len(generated_sample) < min_sample_size:
                    biggest_decrease = {"idx": None, "decrease": 0, "density": 0}
                    for i, density in enumerate(train_densities):
                        try:
                            percent_decrease = (density - gen_densities[i])/density
                            if percent_decrease > biggest_decrease["decrease"] and percent_decrease != 1 and i not in previous:
                                biggest_decrease["idx"] = i
                                biggest_decrease["decrease"] = percent_decrease
                                biggest_decrease["density"] = density
                        except ZeroDivisionError:
                            continue
                    # Selects a random sample from each grouped weight bin for the generated
                    # databases, scaled based on the number of molecules in the selected bin
                    # and the training densities.
                    num_mols_in_bin = len(grouped_data[biggest_decrease["idx"]])
                    for i, density in enumerate(train_densities):
                        sample_len = int(density/biggest_decrease["density"] * num_mols_in_bin)
                        try:
                            grouped_data[i] = random.sample(grouped_data[i], sample_len)
                        except ValueError:
                            continue
                    
                    # Unpacks the sampled weight indexes into one list and adds that list to
                    # another list.
                    generated_sample = [data_idx for data_list in grouped_data.values() for data_idx in data_list]
                    previous.append(biggest_decrease["idx"])

                print(len(generated_sample))
                generated_samples.append(generated_sample)

    # Plots distribution if plot is True.
    if plot is True:
        for i, arg_set in enumerate(args):
            if sample_data == "Molecular Weights":
                _, data = get_atoms_data(arg_set[1].replace("smiles", "ecomp"), use_valid = use_valid)
            else:
                data, _ = get_atoms_data(arg_set[1].replace("smiles", "ecomp"), use_valid = use_valid)
            if i != 0:
                data = [data[idx] for idx in generated_samples[i - 1]]
            label = arg_set[0]

            plt.hist(
                data,
                bins = bins,
                edgecolor = "black",
                width = bin_width,
                rwidth = bin_width,
                label = label,
                alpha = 0.6,
                density = True
                )

        plt.xlabel("Count")
        plt.ylabel("Density")
        # if len(args) == 8:
        #     plt.title(f"Histogram Showing Distribution of {sample_data} in thiols databases")
        # elif len(args) == 2:
        #     plt.title(f"Histogram Showing Distribution of {sample_data} in qm9 databases")
        # elif len(args) == 5:
        #     plt.title(f"Histogram Showing Distribution of {sample_data} in OE62 databases")
        # else:
        #     plt.title(f"Histogram Showing Distribution of {sample_data} in OE62+THz databases")
        plt.margins(x=0)
        plt.legend()
        plt.tight_layout()
        plt.show()
    
        all_data = []
        for i, sample in enumerate(generated_samples):
            if sample_data == "Molecular Weights":
                _, data = get_atoms_data(args[i+1][1].replace("smiles", "ecomp"), use_valid = use_valid)
            else:
                data, _ = get_atoms_data(args[i+1][1].replace("smiles", "ecomp"), use_valid = use_valid)
            sampled_data = [data[i] for i in sample]
            all_data.append(sampled_data)
    
        if sample_data == "Molecular Weights":
            _, original_data = get_atoms_data(args[0][1].replace("smiles", "ecomp"), use_valid = use_valid)
        else:
            original_data, _ = get_atoms_data(args[0][1].replace("smiles", "ecomp"), use_valid = use_valid)
        all_data.insert(0, original_data)

        for i, arg_set in enumerate(args):
            density = gaussian_kde(all_data[i], bw_method=0.3)
            x = np.linspace(0, max(all_data[i]), 1000)
            label = arg_set[0]
            plt.plot(x, density(x), alpha=0.7, label = label)

        plt.ylabel("Density")
        plt.xlabel(sample_data)
        # if len(args) == 8:
        #     plt.title(f"Kernel Density Estimate Showing Distribution of {sample_data} in thiols databases")
        # elif len(args) == 2:
        #     plt.title(f"Kernel Density Estimate Showing Distribution of {sample_data} in qm9 databases")
        # elif len(args) == 5:
        #     plt.title(f"Kernel Density Estimate Showing Distribution of {sample_data} in OE62 databases")
        # else:
        #     plt.title(f"Kernel Density Estimate Showing Distribution of {sample_data} in OE62+THz databases")
        
        plt.margins(x=0)
        plt.legend()
        plt.tight_layout()
        plt.show()
    
    return generated_samples


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
             4. Databases 
            """
            )


        filtered = input("Use filtered data? (y/n)\n")


        if choice1 == "1": 
            if filtered == 'n':## raw 
                args = (
                ("qm9_smiles.db", "qm9_ecomp.json"),
                ("Gschnet_qm9.db", "Gschnet_qm9_ecomp.json"),
                
                ("EDM_qm9.db","EDM_qm9_ecomp.json"),
                ("GeoLDM_qm9.db", "GeoLDM_qm9_ecomp.json")
                )
            else:
                args = (
               ("qm9_smiles.db", "qm9_ecomp.json"),
               ("Gschnet_qm9_filtered.db", "Gschnet_qm9_filtered_ecomp.json"),
               ("EDM_qm9_filtered.db","EDM_qm9_filtered_ecomp.json"),
               ("GeoLDM_qm9_filtered.db","GeoLDM_qm9_filtered_ecomp.json")
               )##filtered

        elif choice1 =="2":
            if filtered == 'y':
                args = (
                    ("geom_drugs.db","geom_drugs_ecomp.json"),
                    ("Gschnet_drugs_filtered.db","Gschnet_drugs_filtered_ecomp.json"),
                
                    ("EDM_drugs_filtered.db","EDM_drugs_filtered_ecomp.json"),
                    ("GeoLDM_drugs_filtered.db","GeoLDM_drugs_filtered_ecomp.json"),
                    ("JODO_drugs_filtered.db","JODO_drugs_filtered_ecomp.json")
                )
            else:
                args = (
                    ("geom_drugs.db","geom_drugs_ecomp.json"),
                    ("Gschnet_drugs_raw.db","Gschnet_drugs_raw_ecomp.json"),

                    ("EDM_drugs_raw.db","EDM_drugs_raw_ecomp.json"),
                    ("GeoLDM_drugs_raw.db","GeoLDM_drugs_raw_ecomp.json"),
                    ("JODO_drugs_raw.db","JODO_drugs_raw_ecomp.json")
                )


        elif choice1 == "3":
            if filtered == 'n':
                args = (
                    ("OE62_full.db", "OE62_full_ecomp.json"),
                    ("Gschnet_oe62_raw.db", "Gschnet_oe62_raw_ecomp.json"),
                    ("EDM_oe62_raw.db", "EDM_oe62_raw_ecomp.json")

                    )
                
            else:
                args = (
                    ("OE62_full.db", "OE62_full_ecomp.json"),
                    ("Gschnet_oe62_filtered.db", "Gschnet_oe62_filtered_ecomp.json"),
                    ("EDM_oe62_filtered.db", "EDM_oe62_filtered_ecomp.json")
                )



        elif choice1 == "4":
            args = (
               ("OE62_full.db", "OE62_full_ecomp.json"),
               ("qm9_filtered.db", "qm9_filtered_ecomp.json"),
                ("geom_drugs.db","geom_drugs_ecomp.json")
            )
            
        

       
        


        

                


       
        else:
            quit()

        # The user can choose to plot data using subsets of the generated
        # databases, sampled such that the molecular weight distribution
        # matches that of the training data.
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
                if choice1 == "12":
                    split_data = input("Split database into OE62 and THz (y/n)?\n")
                    if split_data == "y":
                        split = True
                original_args = [arg_set for arg_set in args]
                generated_samples = sample_generated(args, sample_data = data, bin_width = bin_width, split = split)
                args = original_args
            else:
                generated_samples = None

        # If the user chooses, they can write or rewrite the files containing
        # the element composition data for each database.
        choice2 = input("Generate new element composition data (y/n)?\n")
        if choice2 == "y":
            if choice1 in {"1","2","3"}:
                for arg_set in args:
                    generate_data(arg_set[0], arg_set[1])
            else:
                generate_data(args[0], args[1])

        # The user can choose to either create a histogram or get the average
        # composition of each molecule, as many times as they choose.
        while True:
            choice3 = input(
                """What would you like to do? 
                (enter the corresponding number to choose, or anything else to quit)
                1. Get average molecular composition of molecules by percentage.
                2. Plot distribution of an element across all molecules. 
                3. Plot distribution of number of atoms and molecular weights.
                """
                )
            
            if choice3 not in ("1", "2", "3"):
                break

            '''getting rid of the validity code as moleucles are already filtered based on validity'''
            # choice4 = input(
            #     """Which mode would you like to use?
            #     (enter the corresponding number, or anything else to quit)
            #     1. Use all molecules
            #     2. Use only valid molecules
            #     3. Compare all molecules and valid molecules
            #     """
            #     )
            
            # if choice4 == "1":
            #     mode = "all"
            # elif choice4 == "2":
            #     mode = "valid"
            # elif choice4 == "3":
            #     mode = "compare"
            # else:
            #     break

            mode = "all"

            if choice3 == "1":
                if choice1 in {"1","2","3"}:
                    if not split and choice1 == "12":
                        split_data = input("Split database into OE62 and THz (y/n)?\n")
                        if split_data == "y":
                            split = True
                    x, num_av_counts = create_bar_all(args, mode, generated_samples, split = split)
                    # np.save(f"/root/MChem_DGMs/analysis/Plots/QM9/arrays/y_array_{}",num_av_counts)
                    # np.save(f"/root/MChem_DGMs/analysis/Plots/QM9/arrays/x_array_{}",x)
                else:
                    create_bar_comp(args[0], args[1], mode)

            elif choice3 == "2":
                choice5 = input(
                    """Plot distribution with bins or splines? 
                    (enter corresponding number, or enter anything else to cancel)
                    1. Bins
                    2. Splines
                    3. Both
                    """
                    )
                if choice1 in {"1","2","3"}:
                    if choice5 == "1":
                        create_histogram_all(args, mode, generated_samples)
                    elif choice5 == "2":
                        fit_spline_all(args, mode, generated_samples)
                    elif choice5 == "3":
                        plot_kde_hist_all(args, mode, generated_samples)

                else:
                    if choice5 == "1":
                        create_histogram(args[0], args[1], mode)
                    elif choice5 == "2":
                        fit_spline(args[0], args[1], mode)
                    elif choice5 == "3":
                        plot_hist_kde(args[0], args[1], mode)

            elif choice3 == "3":
                if choice1 in {"1","2","3"}:
                    if not split and choice1 == "12":
                        split_data = input("Split database into OE62 and THz (y/n)?\n")
                        if split_data == "y":
                            split = True
                    plot_atoms_data_all(args, mode, generated_samples, split = split)
                else:
                    plot_atoms_data(args[0], args[1], mode)


if __name__ == "__main__":
    main()

