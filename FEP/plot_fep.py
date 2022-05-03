#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 22 09:05:03 2022

@author: martincup
"""


### load as few packages as necessary
import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

### set directories
dir_work = os.getcwd()
dir_data = dir_work + '/data'
dir_results = dir_work + '/results'

### create function for transforming .xvg data
def parse_XVG(l_filenames=None):

### convert to list if not already
    if not isinstance(l_filenames, list):
        l_filenames = [l_filenames]
        
### instantiate frame for appending
    l_agg = []
    l_labels = []

### raise exception if no file provided
    for filename in l_filenames:
    
### make sure it's an .xvg
        if filename.split('.')[-1].lower() != 'xvg':
            raise Exception('\nFile provided not in .xvg format! Exiting...\n')
            
### otherwise, instantiate array for data        
        else:
            data = []
    
### remove file extension for easier handling
            filename = filename.split('.')[0]
    
### load file as text
            with open(f'{filename}.xvg') as f:
                for line in f:
                    
### get names for plot attributes
                    if ('#' in line.split()[0]) or ('@' in line.split()[0]):
                        if 'xaxis' in line:
                            x = line.split('\"')[-2]
                        elif 'yaxis' in line:
                            y = line.split('\"')[-2]
                        elif ('title' in line) & ('subtitle' not in line):
                            label = line.split('\"')[-2]
                            l_labels.append(label)
                            
### append data to two-dimensional array
                    else:
                        data.append(line.split()[:2])
                        
### transform to frame (try to force to numeric)
                try:
                    df = pd.DataFrame(data, columns=[x,y], dtype=float)
                except:
#                     df = pd.DataFrame(data, columns=[x,y])
                    df = pd.DataFrame(data, dtype=float)

### append frame to list (it's faster this way, trust me)
        df['Legend'] = filename.split('.')[0].split('/')[-1]
        l_agg.append(df)
        
### create aggregate frame
    df_agg = pd.concat(l_agg, ignore_index=True, sort=False)

### print  & save for fidelity
    df_agg.to_csv(f'{filename}_processed.csv', index=False)
        
    return(df_agg, l_labels)

    
### FOR CODE DEVELOPMENT
# l_inputs = ['parse_xvg.py', '-p', 'bar', '-f', 'FEP/barint.xvg']

### get & parse user inputs
l_inputs = sys.argv[1:]

if '-p' not in l_inputs:
    raise Exception('Warning! No plot type detected. Exiting...\n')
idx_plottype = l_inputs.index('-p') + 1
plottype = l_inputs[idx_plottype]

if '-f' not in l_inputs:
    raise Exception('Warning! No filepaths detected. Exiting...\n')
idx_filenames = l_inputs.index('-f') + 1
l_filenames = l_inputs[idx_filenames:]


### parse .xvg data
df, labels = parse_XVG(l_filenames=l_filenames)

### change units & types for FEP data
df[df.columns[0]] = df[df.columns[0]].astype(int) 
df[df.columns[1]] = df[df.columns[1]] * 2.478 # kT --> kJ/mol

### make simple line plot automatically
fig, ax = plt.subplots(1, 1, figsize=(8,5))

### plot single set of data
if df['Legend'].nunique() == 1:
    if plottype.lower() == 'line':
        sns.lineplot(x=df.columns[0], y=df.columns[1], data=df, ax=ax)

    if plottype.lower() == 'bar':
        sns.barplot(x=df.columns[0], y=df.columns[1], data=df, color='white', edgecolor='black', lw=2.5, ax=ax)

    title = labels[0]
    plotname = l_filenames[0].split('.')[0].split('/')[-1] + '_processed'
        
### plot comparison of data (must be of same type!)
else:
    if plottype.lower() == 'line':
        sns.lineplot(x=df.columns[0], y=df.columns[1], hue='Legend', data=df, ax=ax)
        
    if plottype.lower() == 'bar':
        sns.barplot(x=df.columns[0], y=df.columns[1], hue='Legend', data=df, edgecolor='black', lw=2.5, ax=ax)

    title = np.unique(labels)[0]
    plotname = l_filenames[0].split('.')[0].split('/')[-1] + '_comparison'

### make some adjustments to visualization based on data type
if labels[0] == 'Free energy differences':
    xlabel = '$\lambda$'
    ylabel = '$\Delta$G (kJ/mol)'
elif labels[0] == 'Free energy integral':
    xlabel = '$\lambda$'
    ylabel = '$\Delta$G (kJ/mol)'
else:
    xlabel = df.columns[0]
    ylabel = df.columns[1]
    
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)

# ax.set_xlim([df[df.columns[0]].min(), df[df.columns[0]].max()])
ax.set_title(title)

### save to FEP directory
dir_fep = l_filenames[0].split('/')[0]
plt.savefig(f'{dir_fep}/{plotname}.png', bbox_inches='tight', dpi=300)


















