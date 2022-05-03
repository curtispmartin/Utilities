# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

### import modules
import sys

### get protein & ligand(s) structure files
# f_protein = '5w8l_processed.gro'
# l_ligands_ext = ['L0A.gro', 'C0A.gro']

### get protein & ligand(s) structure files
l_inputs = sys.argv
f_protein = sys.argv[1]
if len(l_inputs) - 2 == 1:
    l_ligands_ext = list(sys.argv[2])
elif len(l_inputs) - 2 > 1:
    l_ligands_ext = l_inputs[2:]
else:
    raise Exception('\nNot enough arguments provided! Exiting...\n')

### open protein file
with open(f_protein, 'r') as protein:
    l_prot = protein.read().strip()
    l_prot = l_prot.split('\n')

### loop through ligands
for f_ligand in l_ligands_ext:

### open ligand files
    with open(f_ligand, 'r') as ligand:
        l_lig = ligand.read().strip()
        l_lig = l_lig.split('\n')
    
### insert ligand data into protein data 
    count = int(l_prot[1].strip()) + int(len(l_lig[2:-1]))
    l_prot[1] = f' {count}' # update atom count
    for line in l_lig[2:-1]:
        l_prot.insert(-1, line)

### write to new file
with open('complex.gro', 'w') as plcomplex:
    for line in l_prot:
        plcomplex.write(str(line) + '\n')