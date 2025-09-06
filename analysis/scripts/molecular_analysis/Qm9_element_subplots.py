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

### data goes: ecount, bins 
QM9_C_data= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/data_C_QM9.npy", allow_pickle=True)
Gschnet_C_data= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/data_C_Gschnet.npy", allow_pickle=True)
EDM_C_data= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/data_C_EDM.npy", allow_pickle=True)
GeoLDM_C_data= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/data_C_GeoLDM.npy", allow_pickle=True)

QM9_O_data= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/data_O_QM9.npy", allow_pickle=True)
Gschnet_O_data= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/data_O_Gschnet.npy", allow_pickle=True)
EDM_O_data= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/data_O_EDM.npy", allow_pickle=True)
GeoLDM_O_data= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/data_O_GeoLDM.npy", allow_pickle=True)

QM9_N_data= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/data_N_QM9.npy", allow_pickle=True)
Gschnet_N_data= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/data_N_Gschnet.npy", allow_pickle=True)
EDM_N_data= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/data_N_EDM.npy", allow_pickle=True)
GeoLDM_N_data= np.load("/root/MChem_DGMs/analysis/Plots/QM9/arrays/data_N_GeoLDM.npy", allow_pickle=True)


C_data = np.array([QM9_C_data, Gschnet_C_data, EDM_C_data, GeoLDM_C_data], dtype=object)
O_data = np.array([QM9_O_data, Gschnet_O_data, EDM_O_data, GeoLDM_O_data], dtype=object)
N_data = np.array([QM9_N_data, Gschnet_N_data, EDM_N_data, GeoLDM_N_data], dtype=object)

fig, axs = plt.subplots(3, 1, figsize=(10, 15), sharey=True, sharex=True)


for i in C_data:

    axs[0].hist( x=i[0],bins=i[1], width=0.2 ,edgecolor='black',density=True, rwidth=0.2)
    axs[0].legend(['QM9', 'Gschnet', 'EDM', 'GeoLDM'], loc='upper left', fontsize=size)


for i in O_data:
    axs[1].hist( x=i[0],bins=i[1], width=0.2,edgecolor='black',density=True, rwidth=0.2)

for i in N_data:
    axs[2].hist( x=i[0],bins=i[1], width=0.2,edgecolor='black',density=True, rwidth=0.2)


#plt.show()



axs[0].set_title('Carbon', fontsize=size)
axs[1].set_title('Oxygen', fontsize=size)
axs[2].set_title('Nitrogen', fontsize=size)
axs[0].set_ylabel('Density', fontsize=size)
axs[1].set_ylabel('Density', fontsize=size)
axs[2].set_ylabel('Density', fontsize=size)
axs[2].set_xlabel('Count', fontsize=30)

plt.savefig("/root/MChem_DGMs/analysis/Plots/QM9/arrays/subplots_QM9.pdf", bbox_inches='tight',format='pdf', dpi=300)