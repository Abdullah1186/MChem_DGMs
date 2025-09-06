from math import e
from matplotlib import legend
import matplotlib.pyplot as plt
import matplotlib as mpl
from networkx import density
import numpy as np


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
       #  "axes.prop_cycle": plt.cycler("color", plt.cm.RdYlGn(np.linspace(0,1,8,endpoint=True)))
         'axes.prop_cycle': mpl.cycler(color=["k", "r", "b", "g","y","m"])
                                       }
plt.rcParams.update(params)

nos_of_instances_Alkene= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/nos_of_instances_Alkene.npy", allow_pickle=True) 
proportion_list_Alkene= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/proportion_list_Alkene.npy", allow_pickle=True)
shifted_xs_Alkene= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/shifted_xs_Alkene.npy", allow_pickle=True)

nos_of_instances_OH= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/nos_of_instances_OH.npy", allow_pickle=True) 
proportion_list_OH= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/proportion_list_OH.npy", allow_pickle=True)
shifted_xs_OH= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/shifted_xs_OH.npy", allow_pickle=True)



fig, axs = plt.subplots(2, 1, figsize=(10, 15), sharey=True,sharex=True)



for i in range(4):
        if len(proportion_list_Alkene[i]) < len(nos_of_instances_Alkene):
            new_proportions = np.zeros(len(nos_of_instances_Alkene), dtype = float)
            new_proportions[:len(proportion_list_Alkene[i])] = proportion_list_Alkene[i]
            proportion_list_Alkene[i] = new_proportions

        axs[0].bar(
            shifted_xs_Alkene[i],
            proportion_list_Alkene[i],
            width=0.1,
            edgecolor = "black"
            )
        
        axs[0].legend(['QM9', 'Gschnet', 'EDM', 'GeoLDM'], loc='upper right', fontsize=size)
        

for i in range(4):
        if len(proportion_list_OH[i]) < len(nos_of_instances_OH):
            new_proportions = np.zeros(len(nos_of_instances_OH), dtype = float)
            new_proportions[:len(proportion_list_OH[i])] = proportion_list_OH[i]
            proportion_list_OH[i] = new_proportions

        axs[1].bar(
            shifted_xs_OH[i],
            proportion_list_OH[i],
            width=0.1,
            edgecolor = "black"
            )

axs[0].set_title("Alkene", fontsize=size)
axs[1].set_title("Hydroxy", fontsize=size)
axs[0].set_ylabel("Proportion", fontsize=size)
axs[1].set_ylabel("Proportion", fontsize=size)
axs[1].set_xlabel("Number of instances", fontsize=size)
plt.savefig("/root/MChem_DGMs/analysis/Plots/QM9/arrays/functional_group_analysis.pdf", dpi=300, bbox_inches='tight',format='pdf')