#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 15:25:00 2022

@author: martincup
"""

### import modules
import os
import sys

### get user input...FOR CODE DEVELOPMENT
# l_inputs = ['swap_codes.py', '-i', 'L0A', '-o', '9YA', '-f', 'L0A.mol2', 'L0A.str']

### get user input
l_inputs = sys.argv

### give exception if options not properly defined
if '-f' not in l_inputs:
    raise Exception('\Warning! No filepaths detected. Exiting...\n')
if '-i' not in l_inputs:
    raise Exception('\Warning! No input code detected. Exiting...\n')
if '-o'  not in l_inputs:
    raise Exception('\Warning! No output code detected. Exiting...\n')
    
### parse user input
filenames = l_inputs[l_inputs.index('-f') + 1:] # files for conversion
icode = l_inputs[l_inputs.index('-i') + 1] # input code
ocode = l_inputs[l_inputs.index('-o') + 1] # output code

### loop through files
for filename in filenames:

### start counting # of swaps
    counter = 0

### change filename for replacement
    fsplit = filename.split('.')
    filename_old = f'{fsplit[0]}_old.{fsplit[1]}'
    os.system(f'mv {filename} {filename_old}')
    
### open new file for writing
    with open(filename, 'w') as fnew:
        
### open old file for reading
        with open(filename_old, 'r') as fold:
            for line in fold:
                
### write line w/ replace to new file
                if icode in line:
                    counter += 1
                    fnew.write(line.replace(icode, ocode))
                
### or write original line to new file
                else:
                    fnew.write(line)

### print number of swaps for fidelity
    print(f'Number of {icode} codes swapped for {ocode} in {filename}...\t{counter}')