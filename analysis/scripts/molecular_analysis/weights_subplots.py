from math import e
from matplotlib import legend
import matplotlib.pyplot as plt
import matplotlib as mpl
from networkx import density
import numpy as np
from scipy.stats import gaussian_kde

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

filtered_data=np.load('/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/no_of_atoms_filtered.npy', allow_pickle=True)
unfiltered_data=np.load('/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/no_of_atoms_unfiltered.npy', allow_pickle=True)

fig, ax = plt.subplots(2, 1, figsize=(12, 14),sharey=True, sharex=True)



        
for i in range(5):
        
            density = gaussian_kde(unfiltered_data[0][i], bw_method=0.3)
            x = np.linspace(0, max(unfiltered_data[0][i]), 1000)

            ax[0].plot(x, density(x), alpha=0.7)
            ax[0].legend(['QM9', 'Gschnet', 'EDM', 'GeoLDM', 'JODO'], loc='upper right', fontsize=size)
            ax[0].set_title('Unfiltered Data', fontsize=35)
            ax[0].set_ylabel('Density', fontsize=35)

    

for i in range(5):

            density = gaussian_kde(filtered_data[0][i], bw_method=0.3)
            x = np.linspace(0, max(filtered_data[0][i]), 1000)

            ax[1].plot(x, density(x), alpha=0.7)
            ax[1].set_title('Filtered Data', fontsize=35)
            ax[1].set_ylabel('Density', fontsize=35)
            ax[1].set_xlabel('Count', fontsize=35)

plt.savefig("/root/MChem_DGMs/analysis/Plots/DRUGS/arrays/no_of_atoms_distribution.pdf", dpi=300, bbox_inches='tight',format='pdf')