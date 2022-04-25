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
gmx pdb2gmx -f 5w8l_clean.pdb -o 5w8l_processed.gro -ff charmm36-feb2021 -water tip3p 
gmx pdb2gmx -f 5w8l_monomer.pdb -o 5w8l_monomer.gro -ff charmm36-feb2021 -water tip3p 
python prep_topology.py topol.top 4
# python swap_codes.py -i 'Protein_chain_A     1' -o 'Protein_chain_A     4' -f topol.top

### convert CHARMM to GROMACS formatting
python cgenff_charmm2gmx_py3_nx2.py 9YA 9YA.mol2 9YA.str charmm36-feb2021.ff
python cgenff_charmm2gmx_py3_nx2.py NAI NAI.mol2 NAI.str charmm36-feb2021.ff

### generate GROMACS files for ligands
gmx editconf -f 9YA_ini.pdb -o 9YA.gro
gmx editconf -f NAI_ini.pdb -o NAI.gro

### add ligands to complex
python gen_complex.py 5w8l_processed.gro 9YA.gro NAI.gro

### add ligands to topology
python gen_topology.py topol.top 9YA.gro NAI.gro

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
{ echo "\"System\" & ! a H*"; echo q; } | gmx make_ndx -f 9YA.gro -o index_9YA.ndx
{ echo "\"System\" & ! a H*"; echo q; } | gmx make_ndx -f NAI.gro -o index_NAI.ndx

### generate restraints for ligand(s)
echo 2 | gmx genrestr -f 9YA.gro -n index_9YA.ndx -o posre_9YA.itp -fc 1000 1000 1000
echo 2 | gmx genrestr -f NAI.gro -n index_NAI.ndx -o posre_NAI.itp -fc 1000 1000 1000
python gen_restraints.py topol.top posre_9YA.itp posre_NAI.itp

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

### calculate RMSD over simulation for ligand (9YA)
{ echo 13; echo 13; } | gmx rms -s md.tpr -f md_clean.xtc -o data/md_rmsd_9YA.xvg -tu ns
python parse_xvg.py data/md_rmsd_9YA.xvg

### calculate RMSD over simulation for cofactor (NAI)
{ echo 14; echo 14; } | gmx rms -s md.tpr -f md_clean.xtc -o data/md_rmsd_NAI.xvg -tu ns
python parse_xvg.py data/md_rmsd_NAI.xvg

### calculate distance b/w ligand & cofactor
gmx pairdist -s md.tpr -f md_clean.xtc  -ref 'resname 9YA and name C31 N32 C33 C34 S35' -sel 'resname NAI and name N1N C2N C3N C4N C5N C6N' -o data/dist_9YA_NAI.xvg
python parse_xvg.py data/dist_9YA_NAI.xvg

### calculate distance b/w ligand & residues where halogen bond might exist in P8V
gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resname 9YA and name H06' -sel 'resid 139 and name CB' -o data/dist_9YAH_139.xvg -tu ns
python parse_xvg.py data/dist_9YAH_139.xvg

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
gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resname 9YA and name O38' -sel 'resid 168 and name NH1' -o data/dist_9YA_168-1.xvg -tu ns
gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resname 9YA and name O37' -sel 'resid 168 and name NH2' -o data/dist_9YA_168-2.xvg -tu ns
python parse_xvg.py data/dist_9YA_168-1.xvg data/dist_9YA_168-2.xvg

### count hydrogen bonds b/w protein & ligand
{ echo 1; echo 13; } | gmx hbond -s md.tpr -f md_clean.xtc -num data/md_hbnum.xvg -tu ns
python parse_xvg.py data/md_hbnum.xvg

### calculate solvent-accessible surface area
echo 18 | gmx sasa -s md.tpr -f md_clean.xtc -o data/md_sasa.xvg -b 20 -tu ns
python parse_xvg.py data/md_sasa.xvg

### run quick cluster analysis on backbone
# echo 2 | gmx cluster -s md.tpr -f md_clean.xtc -o data/clust.xpm -g data/clust.log -dist data/clust_dist.xvg -ev data/clust_eig.xvg -clid data/clust_id.xvg
#---------------------------------------------------#















