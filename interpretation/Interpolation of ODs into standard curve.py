#!/usr/bin/env python
# coding: utf-8

# %% Setup


from pprint import pprint
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib
import skimage.io as io
from natsort import natsorted
import seaborn as sns;

sns.set_context("talk")
font = {'size': 10, 'weight': 'normal', 'family': 'arial'}
matplotlib.rc('font', **font)


def antigen2D_to_df1D(xlsx_path, sheet, data_col):
    """
    Convert old 2D output format (per antigen) to 1D dataframe
    :param xlsx_path:
    :param sheet:
    :param data_col:
    :return:
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet, index_col=0)
    df = df.unstack().reset_index(name=data_col)  # unpivot (linearize) the table
    df.rename(columns={'level_1': 'antigen_row', 'level_0': 'antigen_col'}, inplace=True)
    df[['antigen_row', 'antigen_col']] = df[['antigen_row', 'antigen_col']].applymap(int)
    df = df[['antigen_row', 'antigen_col', data_col]]
    df.dropna(inplace=True)
    return df


def well2D_to_df1D(xlsx_path, sheet, data_col):
    """
    Convert new 2D output format (per well) to 1D dataframe
    :param xlsx_path:
    :param sheet:
    :param data_col:
    :return:
    """
    df = pd.read_excel(xlsx_path, sheet_name=sheet, index_col=0)
    df = df.unstack().reset_index(name=data_col)  # unpivot (linearize) the table
    df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
    df['well_id'] = df.row_id + df.col_id.map(str)
    df = df[['well_id', data_col]]
    return df


# %% Set paths

# %% First path
data_folder1 = r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-04-16-04-32-COVID_June4_JBassay_images/Stitched data from multiple pysero outputs CR3022 dilutions/IgG2216'
metadata_path1 = os.path.join(data_folder1, 'pysero_output_data_metadata1.xlsx')
OD_path1 = os.path.join(data_folder1, 'median_ODs.xlsx')
int_path1 = os.path.join(data_folder1, 'median_intensities.xlsx')
bg_path1 = os.path.join(data_folder1, 'median_backgrounds.xlsx')
scienion1_path = os.path.join(data_folder1, '2020-06-04-16-08-27-COVID_June4_JBassay_analysis.xlsx')

# %% Second path
data_folder2 = r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-04-16-04-32-COVID_June4_JBassay_images/Stitched data from multiple pysero outputs CR3022 dilutions/IgG2229'
metadata_path2 = os.path.join(data_folder2, 'pysero_output_data_metadata1.xlsx')
OD_path2 = os.path.join(data_folder2, 'median_ODs.xlsx')
int_path2 = os.path.join(data_folder2, 'median_intensities.xlsx')
bg_path2 = os.path.join(data_folder2, 'median_backgrounds.xlsx')
# scienion_path=os.path.join(data_folder, '2020-05-18-17-59-01-COVID_May18_JVassay_analysis.xlsx')

# %% Third path
data_folder3 = r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-04-16-04-32-COVID_June4_JBassay_images/Stitched data from multiple pysero outputs CR3022 dilutions/Biotin2132'
metadata_path3 = os.path.join(data_folder3, 'pysero_output_data_metadata1.xlsx')
OD_path3 = os.path.join(data_folder3, 'median_ODs.xlsx')
int_path3 = os.path.join(data_folder3, 'median_intensities.xlsx')
bg_path3 = os.path.join(data_folder3, 'median_backgrounds.xlsx')
# scienion_path=os.path.join(data_folder, '2020-05-18-17-59-01-COVID_May18_JVassay_analysis.xlsx')

# %% Fourth path
data_folder4 = r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-04-16-04-32-COVID_June4_JBassay_images/Stitched data from multiple pysero outputs CR3022 dilutions/Biotin2205'
metadata_path4 = os.path.join(data_folder4, 'pysero_output_data_metadata1.xlsx')
OD_path4 = os.path.join(data_folder4, 'median_ODs.xlsx')
int_path4 = os.path.join(data_folder4, 'median_intensities.xlsx')
bg_path4 = os.path.join(data_folder4, 'median_backgrounds.xlsx')
# scienion_path=os.path.join(data_folder, '2020-05-18-17-59-01-COVID_May18_JVassay_analysis.xlsx')

# %% Read antigen and plate info
sheet_names = ['serum ID','serum ID1',
               'serum cat',
               'serum dilution','serum dilution2',
               'serum type',
               'secondary ID',
               'secondary dilution']
plate_info_df = pd.DataFrame()
with pd.ExcelFile(metadata_path1) as metadata_xlsx:
    # get sheet names that are available in metadata
    sheet_names = list(set(metadata_xlsx.sheet_names).intersection(sheet_names))
    for sheet_name in sheet_names:
        sheet_df = pd.read_excel(metadata_path1, sheet_name=sheet_name, index_col=0)
        sheet_df = sheet_df.unstack().reset_index(name=sheet_name)  # unpivot (linearize) the table
        sheet_df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
        if plate_info_df.empty:
            plate_info_df = sheet_df
        else:
            plate_info_df = pd.merge(plate_info_df,
                                     sheet_df,
                                     how='left', on=['row_id', 'col_id'])
plate_info_df['well_id'] = plate_info_df.row_id + plate_info_df.col_id.map(str)
sheet_names.append('well_id')
# convert to number and non-numeric to NaN
plate_info_df['serum dilution'] = \
    plate_info_df['serum dilution'].apply(pd.to_numeric, errors='coerce')
plate_info_df.dropna(inplace=True)
# %%
if np.all(plate_info_df['serum dilution'] >= 1):
    # convert dilution to concentration
    plate_info_df['serum dilution'] = 1 / plate_info_df['serum dilution']
plate_info_df.drop(['row_id', 'col_id'], axis=1, inplace=True)

# %% Read antigen information.
antigen_df = antigen2D_to_df1D(xlsx_path=metadata_path1, sheet='antigen_array', data_col='antigen')

# %% Read optical density from pysero
OD_df = pd.DataFrame()
int_df = pd.DataFrame()
bg_df = pd.DataFrame()
with pd.ExcelFile(OD_path1) as OD_xlsx:
    for _, row in antigen_df.iterrows():
        sheet_name = 'od_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
        OD_1_antiten_df = well2D_to_df1D(xlsx_path=OD_xlsx, sheet=sheet_name, data_col='OD')
        OD_1_antiten_df['antigen_row'] = row['antigen_row']
        OD_1_antiten_df['antigen_col'] = row['antigen_col']
        OD_1_antiten_df['antigen'] = row['antigen']
        OD_df = OD_df.append(OD_1_antiten_df, ignore_index=True)

with pd.ExcelFile(int_path1) as int_xlsx:
    for _, row in antigen_df.iterrows():
        sheet_name = 'int_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
        int_1_antiten_df = well2D_to_df1D(xlsx_path=int_xlsx, sheet=sheet_name, data_col='intensity')
        int_1_antiten_df['antigen_row'] = row['antigen_row']
        int_1_antiten_df['antigen_col'] = row['antigen_col']
        int_df = int_df.append(int_1_antiten_df, ignore_index=True)

with pd.ExcelFile(bg_path1) as bg_xlsx:
    for _, row in antigen_df.iterrows():
        sheet_name = 'bg_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
        bg_1_antiten_df = well2D_to_df1D(xlsx_path=bg_xlsx, sheet=sheet_name, data_col='background')
        bg_1_antiten_df['antigen_row'] = row['antigen_row']
        bg_1_antiten_df['antigen_col'] = row['antigen_col']
        bg_df = bg_df.append(bg_1_antiten_df, ignore_index=True)

# %% merge OD with antigen and plate info.

# Use of filter avoids merge of duplicate columns when the cell is run multiple times.
OD_df = OD_df.filter(items=['antigen_row', 'antigen_col', 'OD', 'well_id'], axis=1)
OD_df = pd.merge(OD_df,
                 antigen_df,
                 how='left', on=['antigen_row', 'antigen_col'])
OD_df = pd.merge(OD_df,
                 plate_info_df,
                 how='right', on=['well_id'])
python_df = pd.merge(OD_df,
                     int_df,
                     how='left', on=['antigen_row', 'antigen_col', 'well_id'])
python_df = pd.merge(python_df,
                     bg_df,
                     how='left', on=['antigen_row', 'antigen_col', 'well_id'])

python_df['pipeline'] = 'python'

# In[171]:
# python_df = python_df.append(scn_df)
python_df.replace([np.inf, -np.inf], np.nan, inplace=True)
python_df.dropna(subset=['OD'], inplace=True)

# %% Do the same for second path
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %% Read antigen and plate info
sheet_names = ['serum ID','serum ID1',
               'serum cat',
               'serum dilution','serum dilution2',
               'serum type',
               'secondary ID',
               'secondary dilution']
plate_info_df2 = pd.DataFrame()
with pd.ExcelFile(metadata_path2) as metadata_xlsx:
    # get sheet names that are available in metadata
    sheet_names = list(set(metadata_xlsx.sheet_names).intersection(sheet_names))
    for sheet_name in sheet_names:
        sheet_df = pd.read_excel(metadata_path2, sheet_name=sheet_name, index_col=0)
        sheet_df = sheet_df.unstack().reset_index(name=sheet_name)  # unpivot (linearize) the table
        sheet_df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
        if plate_info_df2.empty:
            plate_info_df2 = sheet_df
        else:
            plate_info_df2 = pd.merge(plate_info_df2,
                                      sheet_df,
                                      how='left', on=['row_id', 'col_id'])
plate_info_df2['well_id'] = plate_info_df2.row_id + plate_info_df2.col_id.map(str)
sheet_names.append('well_id')
# convert to number and non-numeric to NaN
plate_info_df2['serum dilution'] = \
    plate_info_df2['serum dilution'].apply(pd.to_numeric, errors='coerce')
plate_info_df2.dropna(inplace=True)
# %%
if np.all(plate_info_df2['serum dilution'] >= 1):
    # convert dilution to concentration
    plate_info_df2['serum dilution'] = 1 / plate_info_df2['serum dilution']
plate_info_df2.drop(['row_id', 'col_id'], axis=1, inplace=True)

# %% Read antigen information.
antigen_df2 = antigen2D_to_df1D(xlsx_path=metadata_path2, sheet='antigen_array', data_col='antigen')

# %% Read optical density from pysero
OD_df2 = pd.DataFrame()
int_df2 = pd.DataFrame()
bg_df2 = pd.DataFrame()
with pd.ExcelFile(OD_path2) as OD_xlsx:
    for _, row in antigen_df2.iterrows():
        sheet_name = 'od_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
        OD_1_antiten_df = well2D_to_df1D(xlsx_path=OD_xlsx, sheet=sheet_name, data_col='OD')
        OD_1_antiten_df['antigen_row'] = row['antigen_row']
        OD_1_antiten_df['antigen_col'] = row['antigen_col']
        OD_1_antiten_df['antigen'] = row['antigen']
        OD_df2 = OD_df2.append(OD_1_antiten_df, ignore_index=True)

with pd.ExcelFile(int_path2) as int_xlsx:
    for _, row in antigen_df2.iterrows():
        sheet_name = 'int_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
        int_1_antiten_df = well2D_to_df1D(xlsx_path=int_xlsx, sheet=sheet_name, data_col='intensity')
        int_1_antiten_df['antigen_row'] = row['antigen_row']
        int_1_antiten_df['antigen_col'] = row['antigen_col']
        int_df2 = int_df2.append(int_1_antiten_df, ignore_index=True)

with pd.ExcelFile(bg_path2) as bg_xlsx:
    for _, row in antigen_df2.iterrows():
        sheet_name = 'bg_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
        bg_1_antiten_df = well2D_to_df1D(xlsx_path=bg_xlsx, sheet=sheet_name, data_col='background')
        bg_1_antiten_df['antigen_row'] = row['antigen_row']
        bg_1_antiten_df['antigen_col'] = row['antigen_col']
        bg_df2 = bg_df2.append(bg_1_antiten_df, ignore_index=True)

# %% merge OD with antigen and plate info.


# Use of filter avoids merge of duplicate columns when the cell is run multiple times.
OD_df2 = OD_df2.filter(items=['antigen_row', 'antigen_col', 'OD', 'well_id'], axis=1)
OD_df2 = pd.merge(OD_df2,
                  antigen_df2,
                  how='left', on=['antigen_row', 'antigen_col'])
OD_df2 = pd.merge(OD_df2,
                  plate_info_df2,
                  how='right', on=['well_id'])
python_df2 = pd.merge(OD_df2,
                      int_df2,
                      how='left', on=['antigen_row', 'antigen_col', 'well_id'])
python_df2 = pd.merge(python_df2,
                      bg_df2,
                      how='left', on=['antigen_row', 'antigen_col', 'well_id'])

python_df2['pipeline'] = 'python'
# python_df2.dropna(inplace=True)
# Also update sera type to reflect their identity.
# posseralist=python_df2_fix['Sera ID'].isin(['pos 1','pos 2','pos 3','pos 4'])
# python_df2.loc[posseralist,'type'] = 'Diagnostic'


# In[171]:
# python_df2 = python_df2.append(scn_df)
python_df2.replace([np.inf, -np.inf], np.nan, inplace=True)
python_df2.dropna(subset=['OD'], inplace=True)

#%% Do the same for third path
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%% Read antigen and plate info
sheet_names = ['serum ID','serum ID1',
               'serum cat',
               'serum dilution','serum dilution2',
               'serum type',
               'secondary ID',
               'secondary dilution']
plate_info_df3 = pd.DataFrame()
with pd.ExcelFile(metadata_path3) as metadata_xlsx:
    # get sheet names that are available in metadata
    sheet_names = list(set(metadata_xlsx.sheet_names).intersection(sheet_names))
    for sheet_name in sheet_names:
        sheet_df = pd.read_excel(metadata_path3, sheet_name=sheet_name, index_col=0)
        sheet_df = sheet_df.unstack().reset_index(name=sheet_name)  # unpivot (linearize) the table
        sheet_df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
        if plate_info_df3.empty:
            plate_info_df3 = sheet_df
        else:
            plate_info_df3 = pd.merge(plate_info_df3,
                                      sheet_df,
                                      how='left', on=['row_id', 'col_id'])
plate_info_df3['well_id'] = plate_info_df3.row_id + plate_info_df3.col_id.map(str)
sheet_names.append('well_id')
# convert to number and non-numeric to NaN
plate_info_df3['serum dilution'] = \
    plate_info_df3['serum dilution'].apply(pd.to_numeric, errors='coerce')
plate_info_df3.dropna(inplace=True)
# %%
if np.all(plate_info_df3['serum dilution'] >= 1):
    # convert dilution to concentration
    plate_info_df3['serum dilution'] = 1 / plate_info_df3['serum dilution']
plate_info_df3.drop(['row_id', 'col_id'], axis=1, inplace=True)

# %% Read antigen information.
antigen_df3 = antigen2D_to_df1D(xlsx_path=metadata_path3, sheet='antigen_array', data_col='antigen')

# %% Read optical density from pysero
OD_df3 = pd.DataFrame()
int_df3 = pd.DataFrame()
bg_df3 = pd.DataFrame()
with pd.ExcelFile(OD_path3) as OD_xlsx:
    for _, row in antigen_df3.iterrows():
        sheet_name = 'od_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
        OD_1_antiten_df = well2D_to_df1D(xlsx_path=OD_xlsx, sheet=sheet_name, data_col='OD')
        OD_1_antiten_df['antigen_row'] = row['antigen_row']
        OD_1_antiten_df['antigen_col'] = row['antigen_col']
        OD_1_antiten_df['antigen'] = row['antigen']
        OD_df3 = OD_df3.append(OD_1_antiten_df, ignore_index=True)

with pd.ExcelFile(int_path3) as int_xlsx:
    for _, row in antigen_df3.iterrows():
        sheet_name = 'int_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
        int_1_antiten_df = well2D_to_df1D(xlsx_path=int_xlsx, sheet=sheet_name, data_col='intensity')
        int_1_antiten_df['antigen_row'] = row['antigen_row']
        int_1_antiten_df['antigen_col'] = row['antigen_col']
        int_df3 = int_df3.append(int_1_antiten_df, ignore_index=True)

with pd.ExcelFile(bg_path3) as bg_xlsx:
    for _, row in antigen_df3.iterrows():
        sheet_name = 'bg_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
        bg_1_antiten_df = well2D_to_df1D(xlsx_path=bg_xlsx, sheet=sheet_name, data_col='background')
        bg_1_antiten_df['antigen_row'] = row['antigen_row']
        bg_1_antiten_df['antigen_col'] = row['antigen_col']
        bg_df3 = bg_df3.append(bg_1_antiten_df, ignore_index=True)

# %% merge OD with antigen and plate info.

# Use of filter avoids merge of duplicate columns when the cell is run multiple times.
OD_df3 = OD_df3.filter(items=['antigen_row', 'antigen_col', 'OD', 'well_id'], axis=1)
OD_df3 = pd.merge(OD_df3,
                  antigen_df3,
                  how='left', on=['antigen_row', 'antigen_col'])
OD_df3 = pd.merge(OD_df3,
                  plate_info_df3,
                  how='right', on=['well_id'])
python_df3 = pd.merge(OD_df3,
                      int_df3,
                      how='left', on=['antigen_row', 'antigen_col', 'well_id'])
python_df3 = pd.merge(python_df3,
                      bg_df3,
                      how='left', on=['antigen_row', 'antigen_col', 'well_id'])

python_df3['pipeline'] = 'python'
# python_df3.dropna(inplace=True)\

# In[171]:
#python_df3 = python_df3.append(scn_df)
python_df3.replace([np.inf, -np.inf], np.nan, inplace=True)
python_df3.dropna(subset=['OD'], inplace=True)

# %% Do the same for fourth path
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %% Read antigen and plate info
sheet_names = ['serum ID','serum ID1',
               'serum cat',
               'serum dilution','serum dilution2',
               'serum type',
               'secondary ID',
               'secondary dilution']
plate_info_df4 = pd.DataFrame()
with pd.ExcelFile(metadata_path4) as metadata_xlsx:
    # get sheet names that are available in metadata
    sheet_names = list(set(metadata_xlsx.sheet_names).intersection(sheet_names))
    for sheet_name in sheet_names:
        sheet_df = pd.read_excel(metadata_path4, sheet_name=sheet_name, index_col=0)
        sheet_df = sheet_df.unstack().reset_index(name=sheet_name)  # unpivot (linearize) the table
        sheet_df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
        if plate_info_df4.empty:
            plate_info_df4 = sheet_df
        else:
            plate_info_df4 = pd.merge(plate_info_df4,
                                      sheet_df,
                                      how='left', on=['row_id', 'col_id'])
plate_info_df4['well_id'] = plate_info_df4.row_id + plate_info_df4.col_id.map(str)
sheet_names.append('well_id')
# convert to number and non-numeric to NaN
plate_info_df4['serum dilution'] = \
    plate_info_df4['serum dilution'].apply(pd.to_numeric, errors='coerce')
plate_info_df4.dropna(inplace=True)
# %%
if np.all(plate_info_df4['serum dilution'] >= 1):
    # convert dilution to concentration
    plate_info_df4['serum dilution'] = 1 / plate_info_df4['serum dilution']
plate_info_df4.drop(['row_id', 'col_id'], axis=1, inplace=True)

# %% Read antigen information.
antigen_df4 = antigen2D_to_df1D(xlsx_path=metadata_path4, sheet='antigen_array', data_col='antigen')

# %% Read optical density from pysero
OD_df4 = pd.DataFrame()
int_df4 = pd.DataFrame()
bg_df4 = pd.DataFrame()
with pd.ExcelFile(OD_path4) as OD_xlsx:
    for _, row in antigen_df4.iterrows():
        sheet_name = 'od_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
        OD_1_antiten_df = well2D_to_df1D(xlsx_path=OD_xlsx, sheet=sheet_name, data_col='OD')
        OD_1_antiten_df['antigen_row'] = row['antigen_row']
        OD_1_antiten_df['antigen_col'] = row['antigen_col']
        OD_1_antiten_df['antigen'] = row['antigen']
        OD_df4 = OD_df4.append(OD_1_antiten_df, ignore_index=True)

with pd.ExcelFile(int_path4) as int_xlsx:
    for _, row in antigen_df4.iterrows():
        sheet_name = 'int_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
        int_1_antiten_df = well2D_to_df1D(xlsx_path=int_xlsx, sheet=sheet_name, data_col='intensity')
        int_1_antiten_df['antigen_row'] = row['antigen_row']
        int_1_antiten_df['antigen_col'] = row['antigen_col']
        int_df4 = int_df4.append(int_1_antiten_df, ignore_index=True)

with pd.ExcelFile(bg_path4) as bg_xlsx:
    for _, row in antigen_df4.iterrows():
        sheet_name = 'bg_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
        bg_1_antiten_df = well2D_to_df1D(xlsx_path=bg_xlsx, sheet=sheet_name, data_col='background')
        bg_1_antiten_df['antigen_row'] = row['antigen_row']
        bg_1_antiten_df['antigen_col'] = row['antigen_col']
        bg_df4 = bg_df4.append(bg_1_antiten_df, ignore_index=True)

# %% merge OD with antigen and plate info.


# Use of filter avoids merge of duplicate columns when the cell is run multiple times.
OD_df4 = OD_df4.filter(items=['antigen_row', 'antigen_col', 'OD', 'well_id'], axis=1)
OD_df4 = pd.merge(OD_df4,
                  antigen_df4,
                  how='left', on=['antigen_row', 'antigen_col'])
OD_df4 = pd.merge(OD_df4,
                  plate_info_df4,
                  how='right', on=['well_id'])
python_df4 = pd.merge(OD_df4,
                      int_df4,
                      how='left', on=['antigen_row', 'antigen_col', 'well_id'])
python_df4 = pd.merge(python_df4,
                      bg_df4,
                      how='left', on=['antigen_row', 'antigen_col', 'well_id'])

python_df4['pipeline'] = 'python'
# python_df4.dropna(inplace=True)
# Also update sera type to reflect their identity.
# posseralist=python_df4_fix['Sera ID'].isin(['pos 1','pos 2','pos 3','pos 4'])
# python_df4.loc[posseralist,'type'] = 'Diagnostic'


# In[171]:
# python_df4 = python_df4.append(scn_df)
python_df4.replace([np.inf, -np.inf], np.nan, inplace=True)
python_df4.dropna(subset=['OD'], inplace=True)

#%% Remove failed wells for n dataframes

#%% Remove failed wells for all dataframes
failed_wells = ['B6', 'B8', 'B10', 'B12','F4','H2','H4', 'H6', 'H9', 'H12']
for failed_well in failed_wells:
    python_df = python_df[python_df.well_id != failed_well]

wells2keep = ['B6', 'B8', 'B10', 'B12']
python_df2 = python_df2[(python_df2['well_id'].isin(wells2keep))]

failed_wells3 = ['A6']
for failedwell in failed_wells3:
    python_df3 = python_df3[python_df3.well_id == failedwell]

failed_wells4 = ['A6']
for failedwell in failed_wells4:
    python_df4 = python_df4[python_df4.well_id != failedwell]

#%% Concatenate dataframes
frames = [python_df, python_df2, python_df3,python_df4]
stitchedpython_df = pd.concat(frames)
# ## Fit curves to above plots

import scipy.optimize as optimization
import itertools

def fourPL(x, A, B, C, D):
    return ((A-D)/(1.0+((x/C)**(B))) + D)

def fit2df(df, model, pipeline):
    """fit model to x, y data in dataframe.
    Return a dataframe with fit x, y for plotting
    """
    serum1='CR3022'
    #sera = df[(df['serum ID']==serum1)]
    sera = df['serum ID1'].unique()
    antigens = df['antigen'].unique()
    secondaries = df['secondary ID'].unique()

    keys = itertools.product(sera, antigens, secondaries)
    df_fit = pd.DataFrame(columns=df.columns)
    df_params = pd.DataFrame()
    for serum, antigen, secondary in keys:
        print(serum, antigen)
        sec_dilu_df = df[(df['serum ID1']==serum) &
                    (df['antigen']==antigen) &
                    (df['secondary ID'] == secondary) &
                    (df['pipeline']==pipeline)]
        sec_dilutions = sec_dilu_df['secondary dilution'].unique()
        for sec_dilution in sec_dilutions:
            sub_df = sec_dilu_df[(sec_dilu_df['secondary dilution'] == sec_dilution)].reset_index(drop=True)
            df_fit_temp = pd.DataFrame()
            guess = [0, 1, 5e-4, 1]
            xdata = sub_df['serum dilution2'].to_numpy()
            ydata = sub_df['OD'].to_numpy()
            params, params_covariance = optimization.curve_fit(model, xdata, ydata, guess, bounds=(0, np.inf), maxfev=1e5)
            x_input = np.logspace(np.log10(np.min(xdata)), np.log10(np.max(xdata)), 50)
            y_fit = fourPL(x_input, *params)

            df_params_temp = pd.DataFrame()
            df_params_temp['sc'] = [serum] #for standard curve
            df_params_temp['antigen'] = [antigen]
            df_params_temp['A'] = params[0]
            df_params_temp['B'] = params[1]
            df_params_temp['C'] = params[2]
            df_params_temp['D'] = params[3]

            df_fit_temp['serum dilution2'] = x_input
            df_fit_temp['OD'] = y_fit
            df_fit_temp['serum ID1'] = ' '.join([serum, 'fit'])

            sub_df_expand = pd.concat(
                [sub_df.loc[[0], ['antigen', 'serum ID',
                             'serum type', 'serum cat', 'well_id','serum dilution',
                             'secondary ID',
                             'secondary dilution',
                             'pipeline']]] * len(df_fit_temp.index), axis=0).reset_index(drop=True)
            df_fit_temp = pd.concat([df_fit_temp, sub_df_expand], axis=1)

            df_fit = df_fit.append(df_fit_temp)
            df_params = df_params.append(df_params_temp)
    print('4PL fitting finished')
    return df_fit, df_params

python_df_fit, python_df_params= fit2df(stitchedpython_df, fourPL, 'python')
cr3022=['CR3022']
CR3022_df_fit = python_df_fit[(python_df_fit['serum ID1'] == 'Neg pool 1/200 + CR3022 fit')]
CR3022_df_params = python_df_params[(python_df_params['sc'] =='Neg pool 1/200 + CR3022')]

#%% Append dataframe of CR3022 fit params to OD df
interp_df = stitchedpython_df.merge(CR3022_df_params, left_on='antigen', right_on='antigen', how='outer')

#%% interpolate into curve and add a new column to dataframe
# def fourPLinterp(y, A, B, C, D):
# 	return C*(((y-A)/(D-y))**(1/B))

def fourPLinterp(row):
	return (row['C'])*((((row['OD'])-(row['A']))/((row['D'])-(row['OD'])))**(1/(row['B'])))

interp_df['value on curve'] = interp_df.apply(fourPLinterp,axis=1)
print(interp_df.dtypes)
makereal = np.real(interp_df['value on curve'])
makereal = pd.DataFrame(data=makereal, dtype=np.float64)
makereal.columns = ['value on curve']
print(makereal)
interp_df['value on curve'] =makereal
interp_df.astype({'value on curve': 'float64'}).dtypes
print(interp_df['value on curve'])

# %% Multiply ODs by dilution factor
# interp_df.loc[interp_df['serum ID'] != 'CR3022', 'value on curve'] = interp_df.loc[interp_df['serum ID'] != 'CR3022', 'value on curve'] * 200
# interp_df

# %% Plot the interpolated values

fig_path = os.path.join(data_folder1, 'pysero_plots temp')
os.makedirs(fig_path, exist_ok=True)

#sera_list = natsorted(interp_df['serum ID1'].unique())
sera_list = ['Neg pool','Pos Pool','Neg pool 1/200 + CR3022']
dilutionforplot = ['0.005']
markers = 'o'
sec_dilutions = [2e-4]
hue = 'serum cat'
style = 'serum type'
antigens = natsorted(interp_df['antigen'].unique())
serum_df = interp_df[(interp_df['serum ID1'].isin(sera_list)) & (interp_df['serum dilution'].isin(dilutionforplot))]

assert not serum_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'

g=sns.relplot(x="value on curve", y="OD", hue="serum cat", col_order=antigens, col="antigen", kind="scatter", data=serum_df, col_wrap=5)
g = (g.set_axis_labels("concentration","OD").set(xscale="log"))
plt.savefig(os.path.join(fig_path, 'interp_{}_{}_{}.jpg'.format('positive', 'negative','005')), dpi=300, bbox_inches='tight')

#%% plot the ODs and fit of standard curve
fig_path = os.path.join(data_folder1, 'pysero_plots temp')
os.makedirs(fig_path, exist_ok=True)

pipeline = ['python']
sera_list = ['Neg pool','Pos Pool','Neg pool 1/200 + CR3022']
sera_fit_list = ['Neg pool 1/200 + CR3022 fit']
markers = 'o'
sec_dilutions = [2e-4]
hue = 'serum ID1'
style = 'serum type'
antigens = natsorted(stitchedpython_df['antigen'].unique())

serum_df = stitchedpython_df[(stitchedpython_df['serum ID1'].isin(sera_list))]
assert not serum_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
for sec_id in serum_df['secondary ID'].unique():
    sub_df = serum_df[(serum_df['secondary ID'] == sec_id) & (serum_df['serum dilution']== '0.005')]
    #palette = sns.color_palette(n_colors=len(sub_df[hue].unique()))
    palette = sns.color_palette()
    print('plotting...')
    g = sns.lmplot(x="serum dilution2", y="OD", col_order=antigens, hue=hue, hue_order=sera_list, col="antigen", ci='sd', palette=palette, markers=markers, data=sub_df, col_wrap=5, fit_reg=False, x_estimator=np.mean)
    sub_python_df_fit = python_df_fit[(python_df_fit['serum ID1'].isin(sera_fit_list)) & (python_df_fit['secondary ID'] == sec_id)]
    palette = sns.color_palette(n_colors=len(sub_python_df_fit[hue].unique()))
    for antigen, ax in zip(antigens, g.axes.flat):
        df_fit = sub_python_df_fit[(sub_python_df_fit['antigen'] == antigen)]
        sns.lineplot(x="serum dilution2", y="OD", hue=hue, hue_order=sera_fit_list, data=df_fit,
                     style=style, palette=palette,
                     ax=ax, legend=False)
        ax.set(xscale="log")
        # ax.set(ylim=[-0.05, 1.5])
plt.savefig(os.path.join(fig_path, '{}_{}_{}_fit.jpg'.format('CR3022 005', sec_id, sec_dilutions)),dpi=300, bbox_inches='tight')

# Make combined plot

fig_path = os.path.join(data_folder1, 'pysero_plots temp')
os.makedirs(fig_path, exist_ok=True)

pipeline = ['python']
sera_list = ['Neg pool','Pos Pool']
sera_fit_list = ['Neg pool 1/200 + CR3022 fit']
markers = 'o'
sec_dilutions = [2e-4]
hue = 'serum ID1'
style = 'serum type'
antigens = natsorted(stitchedpython_df['antigen'].unique())

serum_df = interp_df[(interp_df['serum ID1'].isin(sera_list)) & (interp_df['serum dilution'].isin(dilutionforplot))]
assert not serum_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'

for sec_id in serum_df['secondary ID'].unique():
    sub_df = serum_df[(serum_df['secondary ID'] == sec_id) & (serum_df['serum dilution']== '0.005')]
    #palette = sns.color_palette(n_colors=len(sub_df[hue].unique()))
    palette = sns.color_palette()
    print('plotting...')
    g = sns.lmplot(x="value on curve", y="OD", col_order=antigens, hue=hue, hue_order=sera_list, col="antigen", ci='sd', palette=palette, markers=markers, data=sub_df, col_wrap=5, fit_reg=False, x_estimator=np.mean)
    sub_python_df_fit = python_df_fit[(python_df_fit['serum ID1'].isin(sera_fit_list)) & (python_df_fit['secondary ID'] == sec_id)]
    palette = sns.color_palette(n_colors=len(sub_python_df_fit[hue].unique()))
    for antigen, ax in zip(antigens, g.axes.flat):
        df_fit = sub_python_df_fit[(sub_python_df_fit['antigen'] == antigen)]

        sns.relplot(x="value on curve", y="OD", hue="serum cat", col_order=antigens, col="antigen", kind="scatter", data=serum_df, col_wrap=5, ax=ax)
        sns.lineplot(x="serum dilution2", y="OD", hue=hue, hue_order=sera_fit_list, data=df_fit,
                     style=style, palette=palette,
                     ax=ax, legend=False)
        ax.set(xscale="log")
        # ax.set(ylim=[-0.05, 1.5])
plt.savefig(os.path.join(fig_path, '{}_{}_{}_fit.jpg'.format('combined 005', sec_id, sec_dilutions)),dpi=300, bbox_inches='tight')

#%% Plot the positive pool and negative pool curves
fig_path = os.path.join(data_folder1, 'pysero_plots temp')
os.makedirs(fig_path, exist_ok=True)
# serum_type = 'Control'
pipeline = 'python'
sera_list = ['Neg pool','Pos Pool']
sera_fit_list = [' '.join([x, 'fit']) for x in sera_list]
# sera_list = ['pos 1', 'pos 2', 'pos 3', 'pos 4', 'neg 1', 'neg 2', 'neg 3', 'neg 4']
# sera_list = ['CR3022']
# sera_list = ['Neg 1']
# sera_list = ['neg ' + str(i) for i in range(11, 16)]
# sera_list = ['positive ' + str(i) for i in range(1, 6)] + ['negative ' + str(i) for i in range(1, 7)]
# markers = ['o'] * 5 + ['x'] * 6
markers = 'o'
# hue = 'secondary dilution'
hue = 'serum ID1'
sec_dilutions = [2e-4]
# sec_dilutions = [5e-5]
# sec_dilutions = python_df['secondary dilution'].unique()
# style = 'secondary ID'
style = 'serum type'
antigens = natsorted(stitchedpython_df['antigen'].unique())
serum_df = stitchedpython_df[(stitchedpython_df['serum ID1'].isin(sera_list))]
assert not serum_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
for sec_id in serum_df['secondary ID'].unique():
    sub_df = serum_df[(serum_df['secondary ID'] == sec_id)]
    palette = sns.color_palette(n_colors=len(sub_df[hue].unique()))
    print('plotting...')
    g = sns.lmplot(x="serum dilution", y="OD", col_order=antigens,
                    hue=hue, hue_order=sera_list, col="antigen", ci='sd', palette=palette, markers=markers,
                     data=sub_df, col_wrap=5, fit_reg=False, x_estimator=np.mean)
    sub_python_df_fit=python_df_fit[(python_df_fit['pipeline']==pipeline) &
                                   python_df_fit['serum ID1'].isin(sera_fit_list) &
                                    (python_df_fit['secondary ID'] == sec_id) &
                                    python_df_fit['secondary dilution'].isin(sec_dilutions)
                                    ]
    palette = sns.color_palette(n_colors=len(sub_python_df_fit[hue].unique()))
    for antigen, ax in zip(antigens, g.axes.flat):
        df_fit = sub_python_df_fit[(sub_python_df_fit['antigen'] == antigen)]
        sns.lineplot(x="serum dilution", y="OD", hue=hue, hue_order=sera_fit_list, data=df_fit,
                     style=style, palette=palette,
                     ax=ax, legend=False)
        ax.set(xscale="log")
        # ax.set(ylim=[-0.05, 1.5])

    plt.savefig(os.path.join(fig_path, '{}_{}_{}_fit.jpg'.format('_'.join(sera_list), sec_id, sec_dilutions)),
                             dpi=300, bbox_inches='tight')

#%% Test
fig_path = os.path.join(data_folder1, 'pysero_plots temp')
os.makedirs(fig_path, exist_ok=True)

pipeline = ['python']
sera_list = ['Neg pool','Pos Pool']
sera_fit_list = ['Neg pool fit','Pos Pool fit']
markers = 'o'
sec_dilutions = [2e-4]
hue = 'serum ID1'
style = 'serum type'
antigens = natsorted(stitchedpython_df['antigen'].unique())

serum_df = stitchedpython_df[(stitchedpython_df['serum ID1'].isin(sera_list))]
assert not serum_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
for sec_id in serum_df['secondary ID'].unique():
    sub_df = serum_df[(serum_df['secondary ID'] == sec_id)]
    palette = sns.color_palette(n_colors=len(sub_df[hue].unique()))
    #palette = sns.color_palette()
    print('plotting...')
    g = sns.lmplot(x="serum dilution2", y="OD", col_order=antigens, hue=hue, hue_order=sera_list, col="antigen", ci='sd', palette=palette, markers=markers, data=sub_df, col_wrap=5, fit_reg=False, x_estimator=np.mean)
    sub_python_df_fit = python_df_fit[(python_df_fit['serum ID1'].isin(sera_fit_list)) & (python_df_fit['secondary ID'] == sec_id)]
    palette = sns.color_palette(n_colors=len(sub_python_df_fit[hue].unique()))
    for antigen, ax in zip(antigens, g.axes.flat):
        df_fit = sub_python_df_fit[(sub_python_df_fit['antigen'] == antigen)]
        sns.lineplot(x="serum dilution2", y="OD", hue=hue, hue_order=sera_fit_list, data=df_fit,
                     style=style, palette=palette,
                     ax=ax, legend=False)
        ax.set(xscale="log")
        # ax.set(ylim=[-0.05, 1.5])
plt.savefig(os.path.join(fig_path, '{}_{}_{}_fit.jpg'.format('pos neg', sec_id, sec_dilutions)),dpi=300, bbox_inches='tight')