#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 12:22:56 2022

@author: Curtis P. Martin

Parameters
--------------
alpha: p-value considered "statistically significant". Currently only used for plotting
fcthresh: fold-change considered "signficant". log2() transformation applied later

"""

### import modules
import os
import sys

import pandas as pd
import numpy as np
import scipy.stats as stats

import matplotlib.pyplot as plt
import seaborn as sns


##### MAXIMIZE USER INPUT TO IMPROVE WORKFLOW GENERALIZATION
#----------------------------------------------------------------------------#
### set confidence for determining "statistical significance"
alpha = 0.05

### set threshold for determining "high fold-change"
fcthresh = 1.5

### set flag for labeling points (0 means no label; 1 means label)
label_flag = 0

### define data directory... helps keep things organized!
path_data = f'{os.getcwd()}/Data'

### have the user select in the input data
datafiles = pd.Series(os.listdir(path_data))
filename = datafiles.loc[int(input(f'\n{datafiles}\n\nSelect index of input file for analysis: '))].split('.')[0]

### have the user select which sheet in data
datasheets = pd.Series(pd.ExcelFile(f'{path_data}/{filename}.xlsx').sheet_names)
sheetname = datasheets.loc[int(input(f'\n{datasheets}\n\nSelect index of sheet to use for analysis: '))]

### information for loading data
# filename = 'P917-Appendix1-12-IP-samples-MS-data-Summary' # name of spreadsheet (w/o the extension)
# sheetname = 'Summary' # name of sheet of interest w/i file
#----------------------------------------------------------------------------#


##### LOAD, CLEAN, & SELECT DATA
#----------------------------------------------------------------------------#
### load data
# data = pd.read_excel(io=f'{filename}.xlsx', sheet_name=sheetname)
data = pd.read_excel(io=f'{path_data}/{filename}.xlsx', sheet_name=sheetname)

### define & make output directory... helps keep things organized!
path_outp = f'{os.getcwd()}/Output/{filename}'
if not os.path.exists(path_outp):
    os.mkdir(path_outp)

### print out columns for selection by user
l_psmcols1 = input(f'\n{pd.Series(data.columns)}\n\nSelect PSM data for first sample: ').split()
l_psmcols2 = input(f'\n{pd.Series(data.columns)}\n\nSelect PSM data for second sample: ').split()

### get desired name for output files
oname = input('Provide name for output files: ')

### convert indices to integers
l_psmcols1 = [int(i) for i in l_psmcols1]
l_psmcols2 = [int(i) for i in l_psmcols2]

### convert indices to columns
l_psmcols1 = data.columns[l_psmcols1]
l_psmcols2 = data.columns[l_psmcols2]

### make sure data is numeric


### select data of interest
df1 = data[l_psmcols1].copy()
df2 = data[l_psmcols2].copy()

### fill N/A w/ zeros... assuming N/A means not found
# df1 = df1.fillna(0.0)
# df2 = df2.fillna(0.0)
#----------------------------------------------------------------------------#


##### RUN STATISTICAL TESTS
#----------------------------------------------------------------------------#
### create function for generating normalize spectral abundance factor (nSAF)
def calc_nSAF(df, data, col):

### join amino acid counts to frame for metric generation
    df = df.join(data['# AAs']) # '# AAs' might not be standard, so beware!

### generate new metrics 
    df_ret = pd.DataFrame()
    df_ret[f'{col}_AA'] = df[col] / df['# AAs']
    df_ret[f'{col}_SAF'] = df_ret[f'{col}_AA'] / df_ret[f'{col}_AA'].sum()
    df_ret[f'{col}_nSAF'] = df_ret[f'{col}_SAF'] / df_ret[f'{col}_SAF'].max()
    
    return(df_ret)

### normalize PSM --> nSAF
for col in df1.columns:
    df1 = df1.join(calc_nSAF(df1, data, col))
for col in df2.columns:
    df2 = df2.join(calc_nSAF(df2, data, col))
    
### get nSAF columns for statistical tests & plotting later
l_nsafcols1 = [col for col in df1.columns if 'nSAF' in col]
l_nsafcols2 = [col for col in df2.columns if 'nSAF' in col]

### perform two-sample t-test on PSM data
ttest = stats.ttest_ind(a=df1[l_nsafcols1], b=df2[l_nsafcols2], alternative='two-sided', axis=1) # consider `nan_policy='omit'`
l_tstats = ttest[0]
l_pvals = ttest[1]

### join data for further analysis
df_sub = df1.join(df2)
df_sub['p-value'] = l_pvals
#----------------------------------------------------------------------------#


##### PREPARE DATA FOR PLOTTING & OUTPUT
#----------------------------------------------------------------------------#
### filter data to mitochondrial proteins only... MAY NEED TO MAKE MORE FLEXIBLE, DEPENDING ON DATA CONSISTENCY
f_mito = 'Mitochondrial' 
df_vol = data[data['Mitochondrial location'] == f_mito].copy()

### create new frame for volcano plot
df_vol = df_vol[['Accession','GenSymbol']].join(df_sub)

### calculate log2(FC) = log2 of fold-change
# df_vol['Mean_Grp1'] = df_vol[l_psmcols1].mean(axis=1)
# df_vol['Mean_Grp2'] = df_vol[l_psmcols2].mean(axis=1)
df_vol['Mean_Grp1'] = df_vol[l_nsafcols1].mean(axis=1)
df_vol['Mean_Grp2'] = df_vol[l_nsafcols2].mean(axis=1)
df_vol['FC'] = df_vol['Mean_Grp2'] / df_vol['Mean_Grp1'] # interpreting FC to mean from grp1 --> grp2

### calculate values for volcano plot (log2(FC) & log10(p-value))
df_vol['log2(FC)'] = np.log2(df_vol['FC'])
df_vol['-log10(p)'] = -np.log10(df_vol['p-value'])

### add flag for "statistical significance" based on user-defined alpha
# df_vol.loc[df_vol[df_vol['p-value'] <= alpha].index, 'Significance'] = 'Y'
# df_vol.loc[df_vol[df_vol['p-value'] > alpha].index, 'Significance'] = 'N'

### add flag for "high fold-change" based on user-defined threshold
fcthresh_log2 = np.log2(fcthresh)
df_vol.loc[df_vol[np.abs(df_vol['log2(FC)']) > fcthresh_log2].index, 'HighFC'] = 'Yes'
df_vol.loc[df_vol[np.abs(df_vol['log2(FC)']) <= fcthresh_log2].index, 'HighFC'] = 'No'

### sort data for simpler selection of important genes
df_vol = df_vol.sort_values(['HighFC', '-log10(p)'], ascending=False).reset_index(drop=True)
#----------------------------------------------------------------------------#


##### PLOT & SAVE DATA TO FILE
#----------------------------------------------------------------------------#
### set a few plot parameters
xmin = -3
xmax = 3
ymin = 0
ymax = np.max(df_vol['-log10(p)']) + 0.1

### generate plot parameters based on user-defined parameters above
xmin_alpha = 0.20
xmax_alpha = fcthresh_log2 / (xmax - xmin)
ymin_fc = -np.log10(alpha) / (ymax - ymin)
ymax_fc = 0.95

### make a volcano plot
palette = {'Yes':'Blue', 'No':'Silver'}
fig, ax = plt.subplots(1, 1, figsize=(8,6))

sns.scatterplot(x='log2(FC)', y='-log10(p)', hue='HighFC', data=df_vol, palette=palette, s=20, edgecolor=None, lw=0, ax=ax)
# ax.axhline(y=-np.log10(alpha), xmin=xmin_alpha, xmax=0.5-xmax_alpha, color='C3', ls='--')
# ax.axhline(y=-np.log10(alpha), xmin=0.5+xmax_alpha, xmax=1.0-xmin_alpha, color='C3', ls='--')
# ax.axvline(x=-fcthresh_log2, ymin=ymin_fc, ymax=ymax_fc, color='C3', ls='--')
# ax.axvline(x=fcthresh_log2, ymin=ymin_fc, ymax=ymax_fc, color='C3', ls='--')

if label_flag == 1:
    for idx in range(df_vol.shape[0]):
        plt.text(x=0.05+df_vol.loc[idx, 'log2(FC)'], y=df_vol.loc[idx, '-log10(p)'], s=df_vol.loc[idx, 'GenSymbol'], fontdict=dict(color='k', size=5))

ax.set_xlim([xmin, xmax])
ax.set_xlabel('$log_{2}$(Fold Change)')
ax.set_ylim([ymin, ymax])
ax.set_ylabel('$-log_{10}$(p-value)')
# plt.legend(title=f'$log_{2}$(FC) > {fcthresh_log2:0.2f}', loc='upper right')
plt.legend(title=f'Fold change > {fcthresh}', loc='upper right', ncol=2)
plt.title(f'Volcano Plot ({oname})')
plt.savefig(f'{path_outp}/vp_{oname}.png', bbox_inches='tight', dpi=300)
plt.show()

### write data to file
if label_flag == 0:
    df_vol.to_excel(f'{path_outp}/vp_{oname}.xlsx', index=False)
#----------------------------------------------------------------------------#






