#!/bin/bash

# Author: Curtis P. Martin

#---------------------------------------------------#
##### RUN 40NS MD SIMULATION FOR LDHA TETRAMER (6Q13) W/ SINGLE SUBUNIT (A) BOUND TO LIGAND (P8V) & COFACTOR (NAI)
#---------------------------------------------------#

### batch submission parameters
#SBATCH --job-name LDHA_P8V_1_40ns
#SBATCH --mail-type BEGIN,END
#SBATCH --partition=multinode
#SBATCH --constraint=x2695
#SBATCH --time=24:00:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=28
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks=448
#SBATCH --no-requeue
#SBATCH --error=error.txt
#SBATCH --output=output.txt
#SBATCH --exclusive

### load modules
module load gromacs
module load python

### change to & set up data directories
mkdir data
mkdir results


##### PREPARE STRUCTURES FOR SIMULATION
#---------------------------------------------------#
### write topology for protein; use CHARMM36 force field in directory
gmx pdb2gmx -f 6q13_clean.pdb -o 6q13_processed.gro -ff charmm36-feb2021 -water tip3p 
gmx pdb2gmx -f 6q13_monomer.pdb -o 6q13_monomer.gro -ff charmm36-feb2021 -water tip3p 
python prep_topology.py topol.top 4
# python swap_codes.py -i 'Protein_chain_A     1' -o 'Protein_chain_A     4' -f topol.top

### convert CHARMM to GROMACS formatting
python cgenff_charmm2gmx_py3_nx2.py P8V P8V.mol2 P8V.str charmm36-feb2021.ff
python cgenff_charmm2gmx_py3_nx2.py NAI NAI.mol2 NAI.str charmm36-feb2021.ff

### generate GROMACS files for ligands
gmx editconf -f P8V_ini.pdb -o P8V.gro
gmx editconf -f NAI_ini.pdb -o NAI.gro

### add ligands to complex
python gen_complex.py 6q13_processed.gro P8V.gro NAI.gro

### add ligands to topology
python gen_topology.py topol.top P8V.gro NAI.gro

### create unit cell
gmx editconf -f complex.gro -o newbox.gro -bt cubic -d 2.0

### add solvent
gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

### prepare system & add ions
gmx grompp -f ions.mdp -c solv.gro -p topol.top -o ions.tpr
echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

### prepare system, run energy minimization, & generate data for analysis
gmx grompp -f em.mdp -c solv_ions.gro -p topol.top -o em.tpr
mpirun -np $SLURM_NTASKS `which mdrun_mpi` -deffnm em
{ echo Potential; echo 0; } | gmx energy -f em.edr -o data/em_potential.xvg
python parse_xvg.py data/em_potential.xvg
#---------------------------------------------------#


##### RUN NVT EQUILIBRATION
#---------------------------------------------------#
### create index for ligand(s)
{ echo "\"System\" & ! a H*"; echo q; } | gmx make_ndx -f P8V.gro -o index_P8V.ndx
{ echo "\"System\" & ! a H*"; echo q; } | gmx make_ndx -f NAI.gro -o index_NAI.ndx

### generate restraints for ligand(s)
echo P8V | gmx genrestr -f P8V.gro -n index_P8V.ndx -o posre_P8V.itp -fc 1000 1000 1000
echo NAI | gmx genrestr -f NAI.gro -n index_NAI.ndx -o posre_NAI.itp -fc 1000 1000 1000
python gen_restraints.py topol.top posre_P8V.itp posre_NAI.itp

#### generate protein/non-protein index for proper thermostat coupling... LIKELY TO BE PROBLEMATIC
{ echo "1 | 13 | 14"; echo q; } | gmx make_ndx -f em.gro -o index.ndx 

### prepare system & run NVT equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
mpirun -np $SLURM_NTASKS `which mdrun_mpi` -deffnm nvt
{ echo Temperature; echo 0; } | gmx energy -f nvt.edr -o data/nvt_temperature.xvg
python parse_xvg.py data/nvt_temperature.xvg
#---------------------------------------------------#


##### RUN NPT EQUILIBRATION
#---------------------------------------------------#
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
mpirun -np $SLURM_NTASKS `which mdrun_mpi` -deffnm npt
{ echo Pressure; echo 0; } | gmx energy -f npt.edr -o data/npt_pressure.xvg
python parse_xvg.py data/npt_pressure.xvg
#---------------------------------------------------#


##### RUN MD SIMULATION
#---------------------------------------------------#
### prepare system & execute
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -n index.ndx -o md.tpr
mpirun -np $SLURM_NTASKS `which mdrun_mpi` -deffnm md

### clean up trajectory
# { echo Protein; echo System; } | gmx trjconv -s md.tpr -f md.xtc -o md_clean.xtc -center -pbc mol -ur compact
{ echo Protein; echo System; } | gmx trjconv -s md.tpr -f md.xtc -o md_clean.xtc -center -pbc nojump -ur compact
{ echo Protein; echo System; } | gmx trjconv -s md.tpr -f md_clean.xtc -o md_fit.xtc -fit rot+trans

### extract initial frame
echo System | gmx trjconv -s md.tpr -f md_clean.xtc -o md_start.pdb -dump 0
#---------------------------------------------------#


##### DO SOME PRELIMINARY ANALYSIS (FOR WHICH WE CAN USE WHOLE RUN)
#---------------------------------------------------#
### calculate RMSD over simulation for protein
{ echo Backbone; echo Backbone; } | gmx rms -s md.tpr -f md_clean.xtc -o data/md_rmsd.xvg -tu ns
python parse_xvg.py data/md_rmsd.xvg

### calculate RMSD over simulation for ligand
{ echo 13; echo 13; } | gmx rms -s md.tpr -f md_clean.xtc -o data/md_rmsd_P8V.xvg -tu ns
python parse_xvg.py data/md_rmsd_P8V.xvg

### calculate RMSD over simulation for cofactor (NAI)
{ echo 14; echo 14; } | gmx rms -s md.tpr -f md_clean.xtc -o data/md_rmsd_NAI.xvg -tu ns
python parse_xvg.py data/md_rmsd_NAI.xvg

### calculate distance b/w ligand & cofactor
gmx pairdist -s md.tpr -f md_clean.xtc  -ref 'resname P8V and name C18 N19 C20 C21 S22' -sel 'resname NAI and name N1N C2N C3N C4N C5N C6N' -o data/dist_P8V_NAI.xvg
python parse_xvg.py data/dist_P8V_NAI.xvg

### calculate distances b/w important residues
gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resid 102 and name CA' -sel 'resid 238 and name CD2' -o data/dist_102_238.xvg -tu ns
gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resid 103 and name CA' -sel 'resid 241 and name CA' -o data/dist_103_241.xvg -tu ns
gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resid 105 and name CA' -sel 'resid 137 and name CA' -o data/dist_105_137.xvg -tu ns
gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resid 108 and name CA' -sel 'resid 238 and name CA' -o data/dist_108_238.xvg -tu ns

python parse_xvg.py data/dist_102_238.xvg
python parse_xvg.py data/dist_103_241.xvg
python parse_xvg.py data/dist_105_137.xvg
python parse_xvg.py data/dist_108_238.xvg
python parse_xvg.py data/dist_102_238.xvg data/dist_103_241.xvg data/dist_105_137.xvg data/dist_108_238.xvg

### calculate distances b/w ligand & residues
gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resname P8V and name O25' -sel 'resid 168 and name NH1' -o data/dist_P8V_168-1.xvg -tu ns
gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resname P8V and name O24' -sel 'resid 168 and name NH2' -o data/dist_P8V_168-2.xvg -tu ns
python parse_xvg.py data/dist_P8V_168-1.xvg data/dist_P8V_168-2.xvg

### calculate distance of potential halogen bonds
gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resname P8V and name F35' -sel 'resid 139 and name CB' -o data/dist_P8VF_139.xvg -tu ns
python parse_xvg.py data/dist_P8VF_139.xvg

### count hydrogen bonds b/w protein & ligand
{ echo 1; echo 13; } | gmx hbond -s md.tpr -f md_clean.xtc -num data/md_hbnum.xvg -tu ns
python parse_xvg.py data/md_hbnum.xvg

### calculate solvent-accessible surface area
echo 18 | gmx sasa -s md.tpr -f md_clean.xtc -o data/md_sasa.xvg -b 20 -tu ns
python parse_xvg.py data/md_sasa.xvg

### run quick cluster analysis on backbone
# echo 2 | gmx cluster -s md.tpr -f md_clean.xtc -o data/clust.xpm -g data/clust.log -dist data/clust_dist.xvg -ev data/clust_eig.xvg -clid data/clust_id.xvg
#---------------------------------------------------#















