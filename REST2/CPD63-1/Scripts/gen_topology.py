# -*- coding: utf-8 -*-
"""
Mandatory inputs: 
    - Path to topology file ('topol.top')
    - Path(s) to ligand topology file(s) ('ligand.itp')

Optional inputs:
    - '-d' option indicates ligands include duplicates; 
        all topologies after -d are considered duplicates;
        ligand parameter files ('lig.prm') will not be included in topology
"""

### import modules
import os
import sys

### get protein & ligand(s) structure files
# f_topol = 'topol.top'
# l_inputs = ['gen_topology.py', 'topol.top', 'L0A.itp', 'L0B.itp', 'L0C.itp', 'L0D.itp', 'C0A.itp', 'C0B.itp', 'C0C.itp', 'C0D.itp']
# l_inputs = ['gen_topology.py', 'topol.top', 'L0A.itp', 'C0A.itp', '-d', 'L0B.itp', 'L0C.itp', 'L0D.itp', 'C0B.itp', 'C0C.itp', 'C0D.itp']

### get protein & ligand(s) structure files
f_topol = sys.argv[1]
l_inputs = sys.argv
if len(l_inputs) - 2 == 1:
    l_ligands_ext = list(sys.argv[2])
elif len(l_inputs) - 2 > 1:
    l_ligands_ext = l_inputs[2:]
else:
    raise Exception('\nNot enough arguments provided! Exiting...\n')

### check for duplicates flag
if '-d' in l_ligands_ext:
    idx_dup = l_ligands_ext.index('-d')
    l_ligands_ext.pop(idx_dup)
        
### open topology file
with open(f_topol, 'r') as topol:
    l_topol = topol.read().strip().split('\n')

### get indicies for insertions
for line in l_topol:
    if line == '; Include forcefield parameters':
        idx_ffparams = l_topol.index(line)
    if line == '#ifdef POSRES':
        idx_posre = l_topol.index(line)

### create new list w/ truncated ligand codes
l_ligands = []
for ligand in l_ligands_ext:
    ligand = ligand.split('.')[0].upper()
    l_ligands.append(ligand)

### add ligand to molecules list @ end of file... work backwards to maintain indicies
for ligand in l_ligands:
    l_topol.insert(len(l_topol), f'{ligand}     1')

### insert ligand topology
s_topheader = '\n; Include ligand topology'
idx_counter = idx_posre + 3
l_topol.insert(idx_counter, s_topheader)
idx_counter += 1
for ligand in l_ligands:
    s_ligtop = f'#include \"{ligand}.itp\"'
    l_topol.insert(idx_counter, s_ligtop)
    idx_counter += 1

### insert ligand parameters... do not include duplicates if present
s_parheader = '\n; Include  unique ligand parameters'
idx_counter = idx_ffparams + 2
l_topol.insert(idx_counter, s_parheader)
idx_counter += 1
try:
    for ligand in l_ligands[:idx_dup]:
        s_ligpar = f'#include \"{ligand}.prm\"'
        l_topol.insert(idx_counter, s_ligpar)
        idx_counter += 1    
except:
    for ligand in l_ligands:
        s_ligpar = f'#include \"{ligand}.prm\"'
        l_topol.insert(idx_counter, s_ligpar)
        idx_counter += 1


# sys.exit()


### write topology to new file
os.system('mv topol.top topol_old.top')
with open('topol.top', 'w') as newtopol:
    for line in l_topol:
        newtopol.write(str(line) + '\n')
