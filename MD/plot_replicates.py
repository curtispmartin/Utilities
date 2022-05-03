#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 11:07:38 2022

@author: Curtis P. Martin

This script is meant to aggregate replicate MD experiments in a simplem manner. 

Arguments required:
    -f provides the processed .csv data file of interest (from `parse_xvg.py`)
    -e provides the experiment name (should match the directory)
    --omit enables the user to omit replicates from analysis
"""

### import modules
import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

### parse user inputs... FOR CODE DEVELOPMENT ONLY
# l_inputs = ['replicates.py', '-f', 'md_rmsd_NAI_processed.csv', '-e', '01-CPD63', '--omit', '01', '02']

### get user input
l_inputs = sys.argv

### parse user input
filename = l_inputs[l_inputs.index('-f') + 1]
expname = l_inputs[l_inputs.index('-e') + 1]
oname = filename.split('.')[0] + f'_{expname}' 

### give exception if options not properly defined
if '-f' not in l_inputs:
    raise Exception('Warning! No filepath detected. Exiting...\n')
if '-e' not in l_inputs:
    raise Exception('Warning! No experiment name detected. Exiting...\n')

### get replicate directories
replicates = [dir for dir in os.listdir() if not '.' in dir]
replicates.sort()

### allow user to omit replicates
if '--omit' in l_inputs:
    l_omits = l_inputs[l_inputs.index('--omit') + 1:]
    for omit in l_omits:
        replicates.remove(omit)

### load data files 
l_data = []
for rep in replicates:
    path_to_rep = f'{rep}/{expname}/data/{filename}'
    try:
        df = pd.read_csv(path_to_rep)
    except:
        print(f'\nWARNING: {path_to_rep} not found. Moving on...')
        continue
    
### add replicate flag
    df['Replicate'] = rep

### append data to list
    l_data.append(df)

### stack data into single frame
df = pd.DataFrame(columns=l_data[0].columns)
for data in l_data:
    df = pd.concat([df, data])
df = df.reset_index(drop=True)

### make plot w/ some bounds
fig, ax = plt.subplots(1, 1, figsize=(8,5))

sns.lineplot(x=df.columns[0], y=df.columns[1], data=df, ax=ax)

ax.set_xlim([df[df.columns[0]].min(), df[df.columns[0]].max()])
plt.title(f'Replicates ({expname})')
plt.savefig(f'{oname}.png', bbox_inches='tight', dpi=300)
# plt.show()
