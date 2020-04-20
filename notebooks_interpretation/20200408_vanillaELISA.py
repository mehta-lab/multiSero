#!/usr/bin/env python
# coding: utf-8

"""
Experimental notes:
* RBD is coated at 2ug/ml for each well
* The positive samples are serial diluted down each column
* The secondary antibody are serial diluted across rows
"""


#%% Imports

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib
import skimage.io as io
import seaborn as sns

sns.set_context("talk")
font = {'size': 8, 'weight': 'normal', 'family': 'times'}
matplotlib.rc('font', **font)

#%% Set paths

data_folder = r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_leica/20200408_VanillaELISA_Leica2.5x_2'
data_file = os.path.join(data_folder, '4_20_14_37_24/intensities.xlsx')

#%% Read the information as dataframe.


def elisa_xlsx_to_df(xlsx_path, sheet):
    '''
    Read information from a sheet of the excel file organized by the row and column of the plate, and convert it into a columnar dataframe.
    Parameters
    ----------
    xlsx_path: path to excel file
    sheet: string specifying the sheet.

    Returns
    -------
    df_out: tabular dataframe with two columns: 'well_id' and name of the sheet.

    '''
    df = pd.read_excel(xlsx_path, sheet_name=sheet, index_col=0)
    df = df.unstack().reset_index(name=sheet)  # unpivot (linearize) the table
    df.rename(columns={'level_0': 'well_col', 'level_1': 'well_row'}, inplace=True)
    df['well_col'] = df['well_col'].apply(str)
    df['well_id'] = df['well_row'] + df['well_col']
    df=df.set_index('well_id')
    df_out = df[sheet]
    return df_out

sheets = ['antibody_name', 'antibody_concentration', 'sample', 'sample_concentration', 'intensity', 'od']
elisa_df = pd.DataFrame()
for sheet in sheets:
    df_sheet = elisa_xlsx_to_df(data_file, sheet)
    elisa_df=elisa_df.join(df_sheet,how='right')
elisa_df=elisa_df.rename(columns={'antibody_concentration':'ab_dilution',
                                  'antibody_name':'antibody',
                                  'sample_concentration':'sample_dilution'})

#%% Plot OD vs sample dilution per sample per antibody


fig_path = os.path.join(data_folder, 'plots')
os.makedirs(fig_path, exist_ok=True)

samplelist = ['Sample 1','Sample 2','Negative','Blank']
selectsamples_exps = elisa_df[(elisa_df['sample'].isin(samplelist))]
nColors=len(elisa_df['ab_dilution'].unique())
g = sns.relplot(x="sample_dilution", y="od", hue="ab_dilution", col="sample", row='antibody', ci=None,
                data=selectsamples_exps, estimator=np.nanmedian, kind='line', palette=sns.color_palette("RdBu", n_colors=nColors))
g.set(xscale="log",xlim=[1E-4,0.002])
leg=g._legend
leg.set_bbox_to_anchor((1.1,0.5))
plt.savefig(os.path.join(fig_path,'ODvsSampleDilutions.png'),bbox_inches='tight')
plt.show()

#%% Plot OD vs ab dilution per sample per antibody

nColors=len(elisa_df['sample_dilution'].unique())
# Color the lines by unique values of sample dilution.

g = sns.relplot(x="ab_dilution", y="od", hue="sample_dilution", col="sample", row='antibody', ci=None,
                data=selectsamples_exps, estimator=np.nanmedian, kind='line', palette=sns.color_palette("RdBu", n_colors=nColors))
g.set(xscale="log")
leg=g._legend
leg.set_bbox_to_anchor((1.1,0.5))
plt.savefig(os.path.join(fig_path,'ODvsAbDilutions.png'),bbox_inches='tight')
plt.show()