#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MANDATORY INPUT:
    - number of copies to generate (-n)
   - mdp file to be expanded (-f)
"""

### import modules
import os
import sys
import shutil


##### GET & PARSE USER INPUT
#----------------------------------------------------------------------------#
### get user input...FOR CODE DEVELOPMENT
# l_inputs = ['prep_fep.py', '-n', '5', '-s', 'fep.swarm', '-f', 'em.mdp', 'nvt.mdp', 'npt.mdp', 'md.mdp']

### get user input
l_inputs = sys.argv

### give exception if options not properly defined
if '-f' not in l_inputs:
    raise Exception('Warning! No filepaths detected. Exiting...\n')
if '-s' not in l_inputs:
    raise Exception('Warning! No SWARM template detected. Exiting...\n')
if '-n'  not in l_inputs:
    raise Exception('Warning! Number of lambda states unknown. Exiting...\n')

### get number of lambdas in FEP workflow
n_lambda = int(l_inputs[l_inputs.index('-n') + 1]) 

### get files for expansion
l_mdpfiles = l_inputs[l_inputs.index('-f') + 1:] 

### get batch submission file
swarmfile = l_inputs[l_inputs.index('-s') + 1]
swarmfile_core = swarmfile.split('.')[0]
swarmfile_ext = swarmfile.split('.')[1]

### instantiate list of batchfiles for creating swarm file
l_newbatchfiles = []

### make new directory to store FEP MDP files if not already present
path_setup = 'Setup'
if not os.path.exists(path_setup):
    os.mkdir(path_setup)
#----------------------------------------------------------------------------#


##### WRITE NEW MDP FILES FOR EACH LAMBDA
#----------------------------------------------------------------------------#
### create new .mdp file for each lambda
for n in range(n_lambda):
    
### make new directory for lambda if not already present
    path_setup_lambda = f'{path_setup}/{n}'
    if not os.path.exists(path_setup_lambda):
        os.mkdir(path_setup_lambda)
   
### loop through & copy files
    for mdpfile in l_mdpfiles:
    
### get file information 
        mdpfile_core = mdpfile.split('.')[0]
        mdpfile_ext = mdpfile.split('.')[1]    

### make new directory for this part of simulation if not already present
        path_setup_lambda_mdp = f'{path_setup}/{n}/{mdpfile_core.upper()}'
        if not os.path.exists(path_setup_lambda_mdp):
            os.mkdir(path_setup_lambda_mdp)
    
### open MDP files for reading & writing
        with open(mdpfile, 'r') as mdp_read, open(f'{path_setup_lambda_mdp}/{mdpfile_core}_{n}.{mdpfile_ext}', 'w') as mdp_write:
    
### update `init_lambda_state` parameter 
            for line in mdp_read:
                if ('init_lambda_state' in line) & ('=' in line):
                    init_lambda_state = line.split()[-1]
                    line = line.replace(init_lambda_state, str(n)) 

### write line into new file
                mdp_write.write(line)
#----------------------------------------------------------------------------#


##### COPY FILES WHICH NEED NO UPDATING
#----------------------------------------------------------------------------#
### copy other files
l_hidefiles = ['.DS_Store', 'fep.swarm', 'fep-parallel.swarm', 'fep.sh']
for file in os.listdir():
    if (file not in l_inputs) & (len(file.split('.')) > 1) & (file not in l_hidefiles):
        file_core = file.split('.')[0]
        file_ext = file.split('.')[-1]
        for n in range(n_lambda):
            shutil.copy(file, f'{path_setup}/{n}/{file_core}_{n}.{file_ext}')
#----------------------------------------------------------------------------#


##### SET UP SWARM SUBMISSION FILE
#----------------------------------------------------------------------------#
### write swarm file
with open(f'{path_setup}/fep.swarm', 'w') as fep_write:

### start counter 
    counter = 0

### write new swarm commands for each lambda
    for n in range(n_lambda):
#         with open('template-fep.swarm', 'r') as fep_read:
        with open(swarmfile, 'r') as fep_read:
            for line in fep_read:

### catch swarm commands so they only get written once
                if ('#SWARM' in line) & (counter > 0):
                    continue

### otherwise make modifications
                else:

### update file names
                    line = line.replace('.cpt', f'_{n}.cpt') # .cpt files
                    line = line.replace('.gro', f'_{n}.gro') # .gro files
                    line = line.replace('.mdp', f'_{n}.mdp') # .mdp files
                    line = line.replace('.ndx', f'_{n}.ndx') # .ndx files
                    line = line.replace('.top', f'_{n}.top') # .tpr files
                    line = line.replace('.tpr', f'_{n}.tpr') # .tpr files

### update paths
                    line = line.replace('$(pwd)', f'$(pwd)/{n}')
                
### update commands                
                    line = line.replace('em;', f'em_{n};')
                    line = line.replace('nvt;', f'nvt_{n};')
                    line = line.replace('npt;', f'npt_{n};')
                    line = line.replace('md;', f'md_{n};')
                    
                    fep_write.write(line)

### write new line after each lambda... helps identify different processes                    
        fep_write.write('\n')
        
### add to counter to prevent copying swarm commands multiple times               
        counter += 1
#----------------------------------------------------------------------------#






