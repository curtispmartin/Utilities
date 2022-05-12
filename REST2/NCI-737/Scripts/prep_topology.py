# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import sys

### get protein & ligand structure files
f_topol = sys.argv[1]
n_compounds = sys.argv[2]

### for testing purposes
# f_topol = 'topol.top'
# n_compounds = 4

### pull topology file
with open(f_topol, 'r') as topol:
    l_topol = topol.read().strip().split('\n')

### get indicies for insertions
for line in l_topol:

### get molecule type for insertion
    if '[ moleculetype ]' in line:
        idx_moltype = l_topol.index(line)
        moltype = l_topol[idx_moltype + 2].split()[0]

### get index for insertion        
    if '; Compound' in line:
        idx_cpds = l_topol.index(line)

### update number of compounds
idx_cpds = idx_cpds + 1
l_topol.pop(idx_cpds)
l_topol.insert(idx_cpds, f'{moltype}\t{n_compounds}')

### write to new file
os.system('mv topol.top topol_old.top')
with open('topol.top', 'w') as newtopol:
    for line in l_topol:
        newtopol.write(str(line) + '\n')
