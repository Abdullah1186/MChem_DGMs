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
         'axes.prop_cycle': mpl.cycler(color=["k", "r", "b", "g","m"])
                                       }
plt.rcParams.update(params)



H_data = []
C_data = []
H_data.append(np.load("/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/kde_H_GEOM-DRUGS.npy", allow_pickle=True))
H_data.append(np.load("/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/kde_H_Gschnet_filtered.npy", allow_pickle=True))
H_data.append(np.load("/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/kde_H_EDM_filtered.npy", allow_pickle=True))
H_data.append(np.load("/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/kde_H_GeoLDM_filtered.npy", allow_pickle=True))
H_data.append(np.load("/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/kde_H_JODO_filtered.npy", allow_pickle=True))
C_data.append(np.load("/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/kde_C_GEOM-DRUGS.npy", allow_pickle=True))
C_data.append(np.load("/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/kde_C_Gschnet_filtered.npy", allow_pickle=True))
C_data.append(np.load("/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/kde_C_EDM_filtered.npy", allow_pickle=True))
C_data.append(np.load("/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/kde_C_GeoLDM_filtered.npy", allow_pickle=True))
C_data.append(np.load("/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/kde_C_JODO_filtered.npy", allow_pickle=True))

fig, axs = plt.subplots(2, 1, figsize=(12, 15), sharey=True, sharex=True)



for i in H_data:
    axs[0].plot(i[0],i[1], alpha=0.7)
    axs[0].legend(['DRUGS', 'Gschnet_filtered', 'EDM_filtered', 'GeoLDM_filtered', 'JODO_filtered'], loc='upper right', fontsize=size)
    axs[0].set_xlim(0,40)
    axs[0].set_xticks(np.arange(0, 41, 5))


for i in C_data:
    axs[1].plot(i[0],i[1], alpha=0.7)

#plt.show()



axs[0].set_title('Hydrogen', fontsize=size)
axs[1].set_title('Carbon', fontsize=size)
axs[0].set_ylabel('Density', fontsize=size)
axs[1].set_ylabel('Density', fontsize=size)
axs[1].set_xlabel('Count', fontsize=size)

plt.savefig("/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/subplots_DRUGS.pdf", bbox_inches='tight',format='pdf', dpi=300)