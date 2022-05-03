#!/bin/bash

# Author: Curtis P. Martin

#---------------------------------------------------#
##### RUN 40NS MD SIMULATION FOR LDHA TETRAMER (5W8L) W/ SINGLE SUBUNIT (A) BOUND TO LIGAND (9YA) & COFACTOR (NAI)
#---------------------------------------------------#

### batch submission parameters
#SBATCH --job-name LDHA_9YA_1_40ns
#SBATCH --mail-type BEGIN,END
#SBATCH --partition=multinode
#SBATCH --constraint=x2695
#SBATCH --time=48:00:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=28
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks=224
#SBATCH --no-requeue
#SBATCH --error=error.txt
#SBATCH --output=output.txt
#SBATCH --exclusive

### load modules
module load gromacs
module load python

### change to & set up data directories
mkdir FEP
mv Results/*/MD/md_*.xvg FEP/
cd FEP

##### CALCULATE RELATIVE BINDING FREE ENERGIES USING BENNETT'S ACCEPTANCE RATIO
#---------------------------------------------------#
gmx bar -f md_*.xvg -o -oi
cd ..
python plot_fep.py -p bar -f FEP/bar.xvg
python plot_fep.py -p bar -f FEP/barint.xvg
