#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 12:49:57 2022

@author: martincup
"""

### import packages
import sys

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

import statsmodels.formula.api as sm


##### CALL DATA
#----------------------------------------------------------------------------#
### call data
# file = 'MitoVES/CellTiterGlo 22Rv1 MitoVES 72h 2-7-22'
file = '20220317-1_BCA'
data = pd.read_excel(f'{file}.xlsx', index_col='Unnamed: 0')
#----------------------------------------------------------------------------#


##### SET PARAMETERS
#----------------------------------------------------------------------------#
### set volume of sample in wellS
vol_cell = 1 # in ul... IS THIS NEEDED THOUGH? RETHINK.

### set location of 'Reference' column (i.e. BCA standard)
lum_cell = 11

### set dilution concentrations & volumes
vol_dil = 60 # ul
conc_dil = 1 # ug/ul

### set volume percentage desired for lysis buffer addition
conc_lys = 0.20 # at this point, 20% is the default
#----------------------------------------------------------------------------#


##### GENERATE OUTPUTS
#----------------------------------------------------------------------------#
### build regression model for reference curve
model = sm.ols(formula='Concentration ~ Reference', data=data).fit() # concentration in ug/ml (ml NOT ul)
print(model.summary())


### build function w/ regression model
def regression(model=None, lum=None, vol_cell=None):
    if not model:
        raise Exception('No OLS model. Exiting...')
    else:
        conc = (model.params['Intercept'] + model.params['Reference'] * lum) / vol_cell
    
    return(conc)


### plot regression line for visualization... don't save
fig, ax = plt.subplots(1, 1, figsize=(10,6))

sns.regplot(x='Reference', y='Concentration', data=data, ax=ax)

ax.set_title(f'Regression Plot\n({file})\n')
ax.set_xlabel('Reference Signal [Abs]')
ax.set_ylabel('Concentration [ug/ul]')

plt.show()


### apply regression model to generate concentrations
df_conc = data.apply(lambda l: regression(model=model, lum=l, vol_cell=vol_cell))
df_conc = df_conc.loc[:, df_conc.columns != 'Concentration'].copy()
df_conc = df_conc.rename(columns={'Reference':lum_cell})
df_conc[df_conc < 0] = np.nan # convert negatives to zero... shouldn't be necessary


### plot data just cause...
fig, ax = plt.subplots(1, 1, figsize=(10,6))

sns.heatmap(df_conc, annot=True, cmap='RdBu', center=0, fmt='.1f', square=True)

ax.set_title(f'Concentration [ug/ml]\n({file})\n')
ax.xaxis.tick_top()

# plt.savefig(f'{file}_processed.png', bbox_inches='tight', dpi=300)
plt.show()

### save data to file is desired...
# df_conc.to_excel(f'{file}_processed.xlsx')


### create a function for calcuating volume (C1*V1 = C2*V2)
def calc_vol(conc_sol, conc_dil, vol_dil):
    vol_sol = (vol_dil * conc_dil) / conc_sol
    return(vol_sol)


### apply function to transformed data
df_conc_vol = df_conc.apply(lambda l: calc_vol(conc_sol=l, conc_dil=conc_dil, vol_dil=vol_dil))

### plot normalized volumes for comparison
fig, ax = plt.subplots(1, 1, figsize=(10,6))

sns.heatmap(df_conc_vol, annot=True, cmap='RdBu', center=0, fmt='.1f', square=True)

ax.set_title(f'Normalization Volume [ul]\n({file})\n')
ax.xaxis.tick_top()

plt.savefig(f'{file}_processed.png', bbox_inches='tight', dpi=300)
plt.show()

### save data to file is desired...
# df_conc_vol.to_excel(f'{file}_processed_vol.xlsx')


### calculate volume of loading buffer required
vol_load = vol_dil * conc_lys 
df_conc_vol_lys = vol_dil - df_conc_vol - vol_load

### generate means for single samples
df_vol_means = pd.DataFrame.from_dict({'Loading':vol_load, 'Lysis':df_conc_vol_lys.mean(), 'Sample':df_conc_vol.mean()})
df_vol_means['Total'] = df_vol_means.sum(axis=1)
df_vol_means.index.name = 'Columns'
# df_vol_means = df_vol_means.dropna()
print(f'\n\nVolumes needed:\n\n{df_vol_means}')
df_vol_means.T.to_excel(f'{file}_processed.xlsx')
#----------------------------------------------------------------------------#
