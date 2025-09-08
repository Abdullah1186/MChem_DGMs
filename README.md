# Comparison of Molecular Generation Models: GSchNet, EDM, and GeoLDM

## Models 

This repository contains three generative deep learning models for molecular generation:

- [**G-SchNet**](https://github.com/atomistic-machine-learning/schnetpack-gschnet.git)
- [**EDM (E(n)-equivariant Diffusion Model)**](https://github.com/ehoogeboom/e3_diffusion_for_molecules.git)
- [**GeoLDM (Geometric Latent Diffusion Model)**](https://github.com/MinkaiXu/GeoLDM.git)
- [**JODO (Learning Joint 2D & 3D Diffusion Models for Complete Molecule Generation)**](https://github.com/GRAPH-0/JODO.git)

### [Models](/models/) contains all the trained models for each DGM that were trained in, <span style="color: red">this project not including pretrained models.</span>

- EDM and GeoLDM use 
- Gschnet uses 

<span style="color: red"> Something about zip files </span>



## Datasets in this Work


### G-SchNet
- Trained on: `QM9`, `OE62`, and `GEOM-DRUGS`

### EDM
- Trained on: `QM9`, `OE62`, and `GEOM-DRUGS`. EDM repo only trains on QM9 and GEOM-DRUGS, a model that can be trained on the OE62 dataset can be found [here](https://github.com/Abdullah1186/OE62-Trainable-EDM.git)

### GeoLDM
- Trained on: `OE62` only and was not able to train successfully, so no molecules were generated. You can download the GeoLDM publication's pretrained models for `QM9` and `GEOM-DRUGS` from the following link:  🔗 [Download Pretrained Models](https://drive.google.com/drive/folders/1EQ9koVx-GA98kaKBS8MZ_jJ8g4YhdKsL)



## Analysis



- All analysis can be performed on filtered and non-filtered datasets, filtering is done by taking the largest fragment (if disconnected) and removing invalid and non-unique molecules.  
- All analyses (except saturation analysis) can be sampled—either by weight or atom count—to match the training data.  
- The sampling method is the same as used in the [G-SchNet biases paper](https://pubs.acs.org/doi/10.1021/acs.jcim.5c00665) where in most plots data is sampled based on the atom count of the training dataset.

### Downloading all the training datasets and generated molecule datasets needed for analysis.

This data can be downloaded [here](ttps://moleculardatabases.s3.eu-west-2.amazonaws.com/databases.tar.gz).

This is a quick and easy set of commands to download them into the right directory:
```bash 
cd MChem_DGMs/analysis
wget https://moleculardatabases.s3.eu-west-2.amazonaws.com/databases.tar.gz
tar -xzvf databases.tar.gz
```

### Creating SMILES

Some of the analyses run on SMILES through a json file of smiles related to every ASE dataset. All databases in this project already have SMILES in their ASE database and have a json file of SMILES for analysis, however to add SMILES to your own database, and to generate and add a json file of your SMILES to the Database directory for analysis, use the [SMILES](analysis/scripts/database_manipulation/SMILES.py) script.

> Make sure that your database is in the MChem_DGMs/analysis/Database directory

### Filtering datasets

[Filtering](analysis/scripts/filtering) contains scripts and tools for post-processing and filtering generated molecules. 

- [filter_disconnects.py](analysis/scripts/filtering/filter_disconnects.py) will filter the database by taking the largest fragment from a disconnected structure.

- [schnetpack_filter](/analysis/scripts/filtering/schnetpack_filter.sh) will filter the molecules based on validity, uniqueness and radicals, and will print out a .txt file with metrics. This script calls apon the [sort_db](/analysis/scripts/filtering/sort_db.py) script that turns the database into a dictionary, and the [filter_generated](/analysis/scripts/filtering/filter_generated.py) script that filters the dictionary and turns the dictionary back into an ASE database. Finally the [novelty](/analysis/scripts/filtering/novelty.py) script will tell you how many novel molecules you have normally it is negligible. 





### Structural and Molecular Analysis

All scripts in  [molecular_analysis](/analysis/scripts/molecular_analysis/)

###  - Saturation Analysis  
Plots showing the number of H atoms vs the total number of atoms.



###  - Structural Analysis

Includes:
- **Bond distance distributions**





###  - Elemental Distribution Analysis

Includes:
- **Elemental composition**
- **Elemental distribution on an element-by-element basis**
- **Weight distribution**



###  - Functional Group Analysis

Includes:
- **Functional group distribution**
- **Ring chemistry**





## Plots

> **Important!** - [plots](/plots/) contains plots that are used in the thesis. However, due to the compuational expensiveness many of the DRUGS plots will take a while - consider using a HPC. Furthermore, the analysis scripts provided can plot other graphs that are not used in the thesis. 



## GEOM-DRUGS Repository

The GEOM-DRUGS repository explains how the database works, as it is not an ASE-compatible database.

🔗 [Clone the Repo Here](https://github.com/learningmatter-mit/geom.git)



## Disclaimers and Extra Notes

- Many features of the code do not currently work, especially in the functional group script. Feel free to contribute fixes~!
- Francesco Bartucca is the original creator of the elemental, structural, and functional group analysis code. I have only tweaked it for my own data and added some extra features like sampling based on atom count.
- The EDM model for DRUGS is located on the SULIS HPC, which is currently inaccessible. For information on generating molecules using EDM and the DRUGS dataset, see:  
  🔗 [Maurer Group Docs](https://maurergroup.github.io/MaurerGroupDocs/)

---

Feel free to explore, modify, and extend the models and tools provided here for your own molecular generation research~! 💫🧪
