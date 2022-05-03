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
    
Run on shell script as:
    python comp_replicates.py -f <processed_data.csv> -e <experiment1> -c <experiment2> --omit <replicates>
"""

### import modules
import os
import sys
import pandas as pd
import scipy.stats as stats
# import statsmodels.api as sm
# from statsmodels.formula.api import ols
import matplotlib.pyplot as plt
import seaborn as sns

### parse user inputs... FOR CODE DEVELOPMENT ONLY
# l_inputs = ['replicates.py', '-f', 'md_rmsd_NAI_processed.csv', '-e', '01-CPD63', '-c', '02-NCI737', '--omit', '01', '02']
# l_inputs = ['replicates.py', '-f', 'md_rmsf_processed.csv', '-e', '01-CPD63', '-c', '02-NCI737', '--omit', '01', '02']

### get user input
l_inputs = sys.argv

### parse user input
filename = l_inputs[l_inputs.index('-f') + 1]
expname = l_inputs[l_inputs.index('-e') + 1]
compname = l_inputs[l_inputs.index('-c') + 1]
oname = filename.split('.')[0] + f'_{expname}_{compname}' 

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

### function for loading replicate data
def load_data(filename, replicates, expname):

### create empty list for storing replicate frames
    l_data = []

### loop through replicates
    for rep in replicates:

### try to load data
        try:
            path_to_rep = f'{rep}/{expname}/data/{filename}'
            df = pd.read_csv(path_to_rep)

### print warning if data not found
        except:
            print(f'\nWARNING: {path_to_rep} not found. Moving on...')
            continue
    
### label data according to experiment & replicate
        df['Experiment'] = expname
        df['Replicate'] = rep

### append data to list
        l_data.append(df)
    
    return(l_data)

### load data from first experiment
l_data = load_data(filename=filename, replicates=replicates, expname=expname)

### load data from second experiment
l_data_comp = load_data(filename=filename, replicates=replicates, expname=compname)
    
### append second experiment to first
for comp in l_data_comp:
    l_data.append(comp)

### stack data into single frame
df = pd.DataFrame(columns=l_data[0].columns)
for data in l_data:
    df = pd.concat([df, data])
df = df.reset_index(drop=True)

### make sure data is numeric
df[df.columns[0]] = df[df.columns[0]].astype(float)
df[df.columns[1]] = df[df.columns[1]].astype(float)


### make plot w/ some bounds
fig, ax = plt.subplots(1, 1, figsize=(8,5))

sns.lineplot(x=df.columns[0], y=df.columns[1], hue='Experiment', data=df, alpha=0.75, ax=ax)

ax.set_xlim([df[df.columns[0]].min(), df[df.columns[0]].max()])
plt.legend(loc='lower center', ncol=df['Experiment'].nunique())
plt.title(f'{expname} v. {compname}')
plt.savefig(f'{oname}.png', bbox_inches='tight', dpi=300)
# plt.show()


# sys.exit()


### prepare data for statistical tests
data1 = df[df['Experiment'] == df['Experiment'].unique()[0]][df.columns[1]]
data2 = df[df['Experiment'] == df['Experiment'].unique()[1]][df.columns[1]]

### perform Kolmogorov-Smirnov test
p_ks = stats.ks_2samp(data1=data1, data2=data2, alternative='two-sided')[1]
print(f'\nKolmogorov-Smirnov test p-value:\t{p_ks}')

### perform Mann-Whitney U test
p_mw = stats.mannwhitneyu(x=data1, y=data2, alternative='two-sided')[1]
print(f'Mann-Whitney U test p-value:\t\t{p_mw}\n')

### plot distributions
g = sns.displot(x=df.columns[1], hue='Experiment', data=df, kde=True, height=5, aspect=1.33)

g.set(title=f'({expname} v. {compname})')

g.savefig(f'{oname}-distr.png', bbox_inches='tight', dpi=300)


sys.exit()


### generate means & print to csv
df_proc = df.groupby([df.columns[0], 'Experiment']).mean().reset_index()
df_proc.to_csv(f'{oname}.csv', index=False)


sys.exit()


### some ad hoc analysis... testing difference b/w RMSF for E101(A)...
df1_1564 = df[(df['Atom'] == 1564) & (df['Experiment'] == '01-CPD63')]
df2_1564 = df[(df['Atom'] == 1564) & (df['Experiment'] == '02-NCI737')]
stats.ttest_ind(a=df1_1564['(nm)'], b=df2_1564['(nm)'], alternative='two-sided')

### ...& G100(A)...
df1_1579 = df[(df['Atom'] == 1579) & (df['Experiment'] == '01-CPD63')]
df2_1579 = df[(df['Atom'] == 1579) & (df['Experiment'] == '02-NCI737')]
stats.ttest_ind(a=df1_1579['(nm)'], b=df2_1579['(nm)'], alternative='greater')

### ...& E101(B)...
df1_6788 = df[(df['Atom'] == 6788) & (df['Experiment'] == '01-CPD63')]
df2_6788 = df[(df['Atom'] == 6788) & (df['Experiment'] == '02-NCI737')]
stats.ttest_ind(a=df1_6788['(nm)'], b=df2_6788['(nm)'], alternative='less')

### ...& E101(C)...
df1_12012 = df[(df['Atom'] == 12012) & (df['Experiment'] == '01-CPD63')]
df2_12012 = df[(df['Atom'] == 12012) & (df['Experiment'] == '02-NCI737')]
stats.ttest_ind(a=df1_12012['(nm)'], b=df2_12012['(nm)'], alternative='less')

### ...& E101(D)...
df1_17236 = df[(df['Atom'] == 17236) & (df['Experiment'] == '01-CPD63')]
df2_17236 = df[(df['Atom'] == 17236) & (df['Experiment'] == '02-NCI737')]
stats.ttest_ind(a=df1_17236['(nm)'], b=df2_17236['(nm)'], alternative='less')

### try 3-way ANOVA? 
# data = df.rename(columns={'RMSD (nm)':'RMSD', 'Time (ns)':'Time'})
# model = ols('RMSD ~ C(Experiment, Sum) * C(Time, Sum)', data=data).fit()

# aov_table = sm.stats.anova_lm(model, typ=2)
# aov_table










