# PDBBind-Opt Workflow

[![arxiv](https://img.shields.io/badge/arXiv-2411.01223-blue)](https://arxiv.org/abs/2411.01223)

This repository contains scripts of PDBBind-Opt workflow, which organizes a bunch of open-source softwares to probe and fix structural problems in PDBBind.

![workflow](assets/workflow.svg)

## Code availability

+ `pre_process/`: Scripts to prepare PDBBind and BioLiP dataset (identifying ligands and extract binding affinity data)
+ `workflow/`: Codes for PDBBind-Opt worflow
  - `dimorphite_dl`: Package to assign protonation states. We modified the `site_substructures.smarts` to make the rules easier.
  - `fix_ligand.py`: LigandFixer module
  - `fix_protein.py`: ProteinFixer module
  - `process.py`: Main workflow
  - `rcsb.py`: Functions to query RCSB (i.e. downloading files, query SMILES strings)
  - `gather.py`: Functions to create metadata csv files
  - `fix_polymer.py`: Functions to fix polymer ligands
  - `maual_smiles.json`: Manually corrected reference SMILES
  - `building_blocks.csv`: SMILES of alpha-amino acids and common N/C terminal caps. Used to create reference SMILES for polymers
+ `error_fix/`: Contains some error analysis
+ `figshare/`: Metadata of BioLiP2-Opt and PDBBind-Opt dumped in Figshare repo.

## Dataset availability

PDBBind-Opt and BioLiP2-Opt datasets prepared by PDBBind-Opt workflow can be found in this [Figshare repoistory](https://figshare.com/collections/PDBBind_Optimization_to_Create_a_High-Quality_Protein-Ligand_Binding_Dataset_for_Binding_Affinity_Prediction/7520133/1).

## How to reconstruct PDBBind-Opt and BioLiP-Opt

+ **Step 1**: Download PDBBind index file from their official website. Run `download.sh` in the `pre_process` to download BioLiP2 dataset
+ **Step 2**: Run `pre_process/create_dataset_csv.ipynb` to extract binding affinity and identifying ligands. This will give the three csv files
+ **Step 3**: Go to the `workflow` and use the following command to run the workflow
```bash
mkdir ../raw_data
python procees.py -i ../pre_process/BioLiP_bind_sm.csv -d ../raw_data/biolip2_opt
python procees.py -i ../pre_process/PDBBind_poly.csv -d ../raw_data/pdbbind_opt_poly --poly
python procees.py -i ../pre_process/PDBBind_sm.csv -d ../raw_data/pdbbind_opt_sm
```
This will take about one day on a 256-core CPU. If you have more nodes, considering split the input csv file to several chunks and run them in parallel. When the workflow finish, in the output directory, each PDBID will have a folder and if the workflow succeed on this PDBID, there will be a file named `done.tag` under its folder, otherwise ther will be a file named `err`. 
+ **Step 4**: Run the `gather.py` to create metadata files, for example:
```bash
python gather.py -i ../pre_process/BioLiP_bind_sm.csv -d ../raw_data/biolip2_opt -o ../figshare/biolip2_opt/biolip2_opt.csv
```
## Requirements
After `conda create -n PDBBindOPTenv`, most of packages can be directly installed using `pip install`, such as `pip install gemmi`,`pip install rdkit-pypi`, `pip install openmm`. In my experience (HPC, Linux, Python==3.11.9 environment), some packages are not easily installed using `conda install conda-forge` for new people in this area, and they are **openmmforcefields**, **openff**, **pdbfixer** and **openbabel**. 

I recommend [**mamba**](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) (mam**b**a, not mam**d**a). 

- Install [**Miniforge**](https://github.com/conda-forge/miniforge)
  ```
  # in my case, I install
  wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh"
  bash Miniforge3-Linux-x86_64.sh
  ```
- Navigate to `${HOME}` root, you will see new `miniforge3` folder alongside your `miniconda3` folder. In `${HOME}/miniforge3/etc/profile.d/`, you will see `conda.sh` and `mamba.sh`, `source` them
  ```
  source /${HOME}/miniforge3/etc/profile.d/conda.sh
  source /${HOME}/miniforge3/etc/profile.d/mamba.sh
  ```
- At this moment, if we check `conda env list`,we will see
  ```
  # conda environments:
  #
                         /${HOME}/miniconda3
                         /${HOME}/miniconda3/envs/PDBBindOPTenv
  base                   /${HOME}/miniforge3
  ```
- `conda activate /${HOME}/miniconda3/envs/PDBBindOPTenv`
- `mamba install -c conda-forge openmmforcefields`
- `mamba install -c conda-forge openff-toolkit`
- `mamba install -c conda-forge pdbfixer`
- `mamba install -c conda-forge openbabel`
    
  
