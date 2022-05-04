#!/bin/bash

# Author: Curtis P. Martin

#---------------------------------------------------#
##### RUN 5NS MD SIMULATION FOR LDHA TETRAMER (5W8L) W/ SINGLE SUBUNIT (A) BOUND TO LIGAND (9YA) & COFACTOR (NAI) USING REPLICA EXCHANGE W/ SOLUTE SCALING (REST2)
#---------------------------------------------------#

### batch submission parameters
#SBATCH --job-name TEST-160-REST2
#SBATCH --mail-type END
#SBATCH --partition=multinode
#SBATCH --constraint=x2695
#SBATCH --time=48:00:00
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=28
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks=336
#SBATCH --no-requeue
#SBATCH --error=output.txt
#SBATCH --output=output.txt
#SBATCH --exclusive

### load modules
module load gromacs/2018.3+plumed2.5b
module load python


##### PREPARE STRUCTURES FOR SIMULATION
#---------------------------------------------------#
### count number of directories in working
# nrep=$(find . -type d -printf x -maxdepth 1 | wc -c)
# (($nrep--))

### loop through folders for EM
for path in */;
    do

### change directory
        cd $path
        echo $path

### write topology for protein; use CHARMM36 force field in directory
        gmx pdb2gmx -f Structures/5w8l_clean.pdb -o 5w8l_processed.gro -ff charmm36-feb2021 -water tip3p 
        gmx pdb2gmx -f Structures/5w8l_monomer.pdb -o 5w8l_monomer.gro -ff charmm36-feb2021 -water tip3p 
        python Scripts/prep_topology.py topol.top 4

### convert CHARMM to GROMACS formatting
        python Scripts/cgenff_charmm2gmx_py3_nx2.py 9YA Structures/9YA.mol2 Structures/9YA.str charmm36-feb2021.ff
        python Scripts/cgenff_charmm2gmx_py3_nx2.py NAI Structures/NAI.mol2 Structures/NAI.str charmm36-feb2021.ff

### generate GROMACS files for ligands
        gmx editconf -f 9YA_ini.pdb -o 9YA.gro
        gmx editconf -f NAI_ini.pdb -o NAI.gro

### add ligands to complex
        python Scripts/gen_complex.py 5w8l_processed.gro 9YA.gro NAI.gro

### add ligands to topology
        python Scripts/gen_topology.py topol.top 9YA.gro NAI.gro

### create unit cell
        gmx editconf -f complex.gro -o newbox.gro -bt cubic -d 1.0

### add solvent
        gmx solvate -cp newbox.gro -cs spc216.gro -p topol.top -o solv.gro

### prepare system & add ions
        gmx grompp -f Setup/ions.mdp -c solv.gro -p topol.top -o ions.tpr
        echo SOL | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

### prepare energy minimization (& prepare pre-processed topology for REST2 scaling)
#         gmx grompp -f Setup/em.mdp -c solv_ions.gro -p topol.top -pp processed.top -o em.tpr
	
### mark "hot" atoms for REST2 (i.e. the solute)
#         python Scripts/mark_hottop.py processed.top Protein_chain_A 9YA NAI

### scale pre-processed topology for REST2
#         plumed partial_tempering 1.0 < processed.top > scaled.top 
#         mv scaled.top scaled_bad.top
#         head -n -11 scaled_bad.top > scaled.top # DUMB
		
### regenerate .tpr file using scaled topology
#         gmx grompp -f Setup/em.mdp -c solv_ions.gro -p scaled.top -o em.tpr
        gmx grompp -f Setup/em.mdp -c solv_ions.gro -p topol.top -o em.tpr
        cd ..

    done

### run multi-directory minimization in preparation for REST2
mpirun -np $SLURM_NTASKS `which mdrun_mpi` -deffnm em -multidir 0[012] # AUTOMATE REPLICA DEFINITION
#---------------------------------------------------#


##### RUN NVT EQUILIBRATION
#---------------------------------------------------#
### loop through folders for NVT
for path in */;
    do

### change directory
        cd $path
        echo $path

### create index for ligand(s)
        { echo "\"System\" & ! a H*"; echo q; } | gmx make_ndx -f 9YA.gro -o index_9YA.ndx
        { echo "\"System\" & ! a H*"; echo q; } | gmx make_ndx -f NAI.gro -o index_NAI.ndx

### generate restraints for ligand(s)
        echo 2 | gmx genrestr -f 9YA.gro -n index_9YA.ndx -o posre_9YA.itp -fc 1000 1000 1000
        echo 2 | gmx genrestr -f NAI.gro -n index_NAI.ndx -o posre_NAI.itp -fc 1000 1000 1000
        python Scripts/gen_restraints.py topol.top posre_9YA.itp posre_NAI.itp

#### generate protein/non-protein index for proper thermostat coupling... LIKELY TO BE PROBLEMATIC
        { echo "1 | 13 | 14"; echo q; } | gmx make_ndx -f em.gro -o index.ndx 

### prepare NVT equilibration
#         gmx grompp -f Setup/nvt.mdp -c em.gro -r em.gro -p topol.top -pp processed.top -n index.ndx -o nvt.tpr

### mark "hot" atoms for REST2 (i.e. the solute)
#         python Scripts/mark_hottop.py processed.top Protein_chain_A 9YA NAI

### scale pre-processed topology for REST2
#         plumed partial_tempering 1.0 < processed.top > scaled.top 
#         mv scaled.top scaled_bad.top
#         head -n -11 scaled_bad.top > scaled.top # DUMB

### regenerate .tpr file using scaled topology
#         gmx grompp -f Setup/nvt.mdp -c em.gro -r em.gro -p scaled.top -n index.ndx -o nvt.tpr
        gmx grompp -f Setup/nvt.mdp -c em.gro -r em.gro -p topol.top -n index.ndx -o nvt.tpr
        cd ..

    done

### run multi-directory NVT in preparation for replica exchange
mpirun -np $SLURM_NTASKS `which mdrun_mpi` -deffnm nvt -multidir 0[012]
#---------------------------------------------------#


##### RUN NPT EQUILIBRATION
#---------------------------------------------------#
### loop through folders for NPT
for path in */;
    do

### change directory
        cd $path

### prepare NPT equilibration
        gmx grompp -f Setup/npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr

### mark "hot" atoms for REST2 (i.e. the solute)
#         python Scripts/mark_hottop.py processed.top Protein_chain_A 9YA NAI

### scale pre-processed topology for REST2
#         plumed partial_tempering 1.0 < processed.top > scaled.top 
#         mv scaled.top scaled_bad.top
#         head -n -11 scaled_bad.top > scaled.top # DUMB
		
### regenerate .tpr file using scaled topology
#         gmx grompp -f Setup/npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p scaled.top -n index.ndx -o npt.tpr
        gmx grompp -f Setup/npt.mdp -c nvt.gro -t nvt.cpt -r nvt.gro -p topol.top -n index.ndx -o npt.tpr
        cd ..

    done

### run multi-directory NPT in preparation for REST2
mpirun -np $SLURM_NTASKS `which mdrun_mpi` -deffnm npt -multidir 0[012]
#---------------------------------------------------#


##### RUN REPLICA EXCHANGE W/ SOLUTE SCALING (REST2)
#---------------------------------------------------#
### loop through folders for REST2
for path in */;
    do

### change directory
        cd $path

### prepare REST2 execution
        gmx grompp -f Setup/md.mdp -c npt.gro -t npt.cpt -p topol.top -pp processed.top -n index.ndx -o md.tpr

### mark "hot" atoms for REST2 (i.e. the solute)
        python Scripts/mark_hottop.py processed.top Protein_chain_A 9YA NAI

### scale pre-processed topology for REST2
        plumed partial_tempering 1.0 < processed.top > scaled.top # AUTOMATE SCALING FACTOR
        mv scaled.top scaled_bad.top
        head -n -11 scaled_bad.top > scaled.top # DUMB

### regenerate .tpr file using scaled topology
        gmx grompp -f Setup/md.mdp -c npt.gro -t npt.cpt -p scaled.top -n index.ndx -o md.tpr
        cd ..

    done

### run multi-directory REST2 
mpirun -np $SLURM_NTASKS `which mdrun_mpi` -deffnm md -plumed Setup/plumed.dat -multidir 0[012] -replex 160 -hrex -dlb no -dhdl dhdl.xvg
#---------------------------------------------------#


##### RUN SOME PRELIMINARY DATA ANALYSIS
#---------------------------------------------------#
### loop through folders for analysis
for path in */;
    do
    
### clean up trajectory for analysis
        { echo Protein; echo System; } | gmx trjconv -s md.tpr -f md.xtc -o md_clean.xtc -center -pbc nojump -ur compact
        { echo Protein; echo System; } | gmx trjconv -s md.tpr -f md_clean.xtc -o md_fit.xtc -fit rot+trans

### calculate RMSD over simulation for protein
        { echo Backbone; echo Backbone; } | gmx rms -s md.tpr -f md_clean.xtc -o md_rmsd.xvg -tu ns
        python parse_xvg.py md_rmsd.xvg

### calculate RMSD over simulation for ligand (9YA)
        { echo 13; echo 13; } | gmx rms -s md.tpr -f md_clean.xtc -o md_rmsd_9YA.xvg -tu ns
        python parse_xvg.py md_rmsd_9YA.xvg

### calculate RMSD over simulation for cofactor (NAI)
        { echo 14; echo 14; } | gmx rms -s md.tpr -f md_clean.xtc -o md_rmsd_NAI.xvg -tu ns
        python parse_xvg.py md_rmsd_NAI.xvg

### calculate distance b/w ligand & cofactor
        gmx pairdist -s md.tpr -f md_clean.xtc  -ref 'resname 9YA and name C31 N32 C33 C34 S35' -sel 'resname NAI and name N1N C2N C3N C4N C5N C6N' -o dist_9YA_NAI.xvg
        python parse_xvg.py dist_9YA_NAI.xvg

### calculate distance b/w ligand & residues where halogen bond might exist in P8V
        gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resname 9YA and name H06' -sel 'resid 139 and name CB' -o dist_9YAH_139.xvg -tu ns
        python parse_xvg.py dist_9YAH_139.xvg

### calculate distances b/w important residues
        gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resid 102 and name CA' -sel 'resid 238 and name CD2' -o dist_102_238.xvg -tu ns
        gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resid 103 and name CA' -sel 'resid 241 and name CA' -o dist_103_241.xvg -tu ns
        gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resid 105 and name CA' -sel 'resid 137 and name CA' -o dist_105_137.xvg -tu ns
        gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resid 108 and name CA' -sel 'resid 238 and name CA' -o dist_108_238.xvg -tu ns
        
        python parse_xvg.py dist_102_238.xvg
        python parse_xvg.py dist_103_241.xvg
        python parse_xvg.py dist_105_137.xvg
        python parse_xvg.py dist_108_238.xvg
        python parse_xvg.py dist_102_238.xvg dist_103_241.xvg dist_105_137.xvg dist_108_238.xvg

### calculate distances b/w ligand & residues
        gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resname 9YA and name O38' -sel 'resid 168 and name NH1' -o dist_9YA_168-1.xvg -tu ns
        gmx pairdist -s md.tpr -f md_clean.xtc -ref 'resname 9YA and name O37' -sel 'resid 168 and name NH2' -o dist_9YA_168-2.xvg -tu ns
        python parse_xvg.py dist_9YA_168-1.xvg dist_9YA_168-2.xvg

### count hydrogen bonds b/w protein & ligand
        { echo 1; echo 13; } | gmx hbond -s md.tpr -f md_clean.xtc -num md_hbnum.xvg -tu ns
        python parse_xvg.py md_hbnum.xvg

### calculate solvent-accessible surface area
        echo 18 | gmx sasa -s md.tpr -f md_clean.xtc -o md_sasa.xvg -tu ns
        python parse_xvg.py md_sasa.xvg
        cd ..
        
    done
#---------------------------------------------------#


##### CLEAN UP THE DIRECTORY
#---------------------------------------------------#
### loop through folders for clean up
for path in */;
    do

### move data to directory
        mkdir Data
        mv *.csv Data/
        mv *.log Data/
        mv *.xvg Data/

### move results to directory
        mkdir Results
        mv *.png Results/
        mv *.txt Results/

### move simulation files to directory & zip them up for simpler IO
        mkdir Output
        mv *.cpt Output/
        mv *.edr Output/
        mv *.gro Output/
        mv *.itp Output/
        mv *.ndx Output/
        mv *.prm Output/
        mv *.top Output/
        mv *.trr Output/
        
        tar -czvf output.tar.gz Output/
        rm -r Output/
        cd ..
        
    done
#---------------------------------------------------#

