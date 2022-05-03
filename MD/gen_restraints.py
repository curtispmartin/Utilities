# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

### import modules
import os
import sys

### get protein & ligand(s) structure files
# f_topol = 'topol.top'
# l_inputs = ['gen_restraints.py', 'topol.top', 'L0A.itp', 'L0B.itp', 'L0C.itp', 'L0D.itp', 'C0A.itp', 'C0B.itp', 'C0C.itp', 'C0D.itp']

### get protein & ligand(s) structure files
f_topol = sys.argv[1]
l_inputs = sys.argv
if len(l_inputs) - 2 == 1:
    l_ligands_ext = list(sys.argv[2])
elif len(l_inputs) - 2 > 1:
    l_ligands_ext = l_inputs[2:]
else:
    raise Exception('\nNot enough arguments provided! Exiting...\n')

### create new list w/ truncated ligand codes
l_ligands = []
for ligand in l_ligands_ext:
    ligand = ligand.split('.')[0].split('_')[-1].upper()
    l_ligands.append(ligand)

### open topology file
with open(f_topol, 'r') as topol:
    l_topol = topol.read().strip().split('\n')

### get indicies for insertions
for line in l_topol:
    if line == '; Include ligand topology':
        idx_ligtopol = l_topol.index(line)

### insert ligand restraints
s_ligrestrheader = '\n; Include position restraints\n#ifdef POSRES'
idx_counter = idx_ligtopol + len(l_inputs) - 1 # len(l_inputs) = restraint files + 2 (.py & .top files)
l_topol.insert(idx_counter, s_ligrestrheader)
idx_counter += 1
for ligand in l_ligands:
    s_ligrestr = f'#include \"posre_{ligand}.itp\"'
    l_topol.insert(idx_counter, s_ligrestr)
    idx_counter += 1
s_ligrestrfooter = '#endif'
l_topol.insert(idx_counter, s_ligrestrfooter)  

### write to new file
os.system('mv topol_old.top topol_older.top')
os.system('mv topol.top topol_old.top')
with open('topol.top', 'w') as newtopol:
    for line in l_topol:
        newtopol.write(str(line) + '\n')
