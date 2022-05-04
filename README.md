# Utilities
Some code I've written & attempted to generalize for future usage. 

## Assays
Contains scripts developed to speed up cell assay analyses. 
- `assay_BCA.py` translates intensity data generated by BCA assay (or similar) & converts to concentrations. Also provides dilutions if/when necessary. 
- `assay_BCA_template.py` provides template needed for running `assay_BCA.py` -- Plan experiment accordingly! 

-----

## FEP
Contains scripts written to execute free energy perturbations (FEP) w/ replica exhange (of some variant). WIP. Not even close to being ready. 

-----

## MD
Contains multiple scripts written to enable simple molecular dynamics (MD) simulations using GROMACS. 
- `cgenff_charmm2gmx_py3_nx2.py` converts ligand parameters (from CGenFF) to GROMACS format. Slightly adapted from [the Mackerell lab](mackerell.umaryland.edu/download.php?filename=CHARMM_ff_params_files/cgenff_charmm2gmx.py).
- `comp_replicates.py` compares data generated from two different experiments, both done in replicate. WIP. 
- `gen_complex.py` creates a `complex.gro` from `protein.gro` & `ligand.gro`, which are processed using `pdb2gmx`
- `gen_restraints.py` inserts all ligand restraints into `topol.top`.
- `gen_topology.py` inserts ligand topologies into `topol.top`.
- `parse_xvg.py` generates simple lineplots from resulting .xvg files produced by GROMACS.
- `plot_replicates.py` starts analysis of replicates in single experiment. WIP. 
- `simulation_*.sh` are the shell scripts for submitting simulation to HPC via `sbatch`.
- `swap_codes.py` enables quick swapping of ligand/protein codes (or anything else really) w/i files. Simplifies simulation setup significantly (somewhat surprisingly, given it's so simple). Provide input code using `-i`, output code using `-o`, & all files for which you want to swap the two using `-f`. 
- `*.mdp` scripts provide the parameters for each part of the simulation. 

-----

## MS
Contains scripts developed to analyze mass spectrometry (MS) data in a user-friendly way (as much as possible w/o creating an .exe).
- `analysis_MS.py` creates volcano plots using MS data. Requires multiple user inputs which are requested interactively during execution. 
- `msenv.yml` provides a virtual environment used for this code. Add virtual enviornment using `conda env create --file=msenv.yml`

-----

## REST2
Contains all scripts required to run replica exchange w/ solute scaling (REST2) simulations in GROMACS. NOTE: Plumed extension is required. 
- `mark_hottop.py` marks atoms desired to participate in replica exchange using Plumed syntax. Requires definition of molecules of interest; right now limited to complete molecules. 
- `sim-rest2.sh` implements REST2 in HPC via SLURM. 
All other code pulled (or slightly modified) from **MD** scripts. 

-----

## RFR
Code written to enable simple implementation of a random forest regression model adapted to time-series forecasting. 
- `fxn_RFR.py` enables simple constructions & validation of a time-series random forest model using the `build_model` class. 
- `env_RFR.yml` provides a virtual environment for simple use of code. Just clone & go. 

-----

## SQL
Code written to enable simple pulls of data from SQL servers. 
- `fxn_SQL.py` provides a `pull_SQL()` function which acts as a wrapper for simple SQL pulls. Great for folks who don't know SQL well but need to access server quickly... Will need slight modifications if manual user credentials required. 
- `env_SQL.py` provides a virtual environment for simple use of code. Just clone & go.

-----