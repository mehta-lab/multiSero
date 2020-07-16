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
sheet_names = ['serum ID',
               'serum dilution',
               'serum cat',
               'serum type',
               'secondary ID',
               'secondary dilution']
plate_info_df = pd.DataFrame()
with pd.ExcelFile(metadata_path1) as metadata_xlsx:
    # get sheet names that are available in metadata
    sheet_names = list(set(metadata_xlsx.sheet_names).intersection(sheet_names))
    for sheet_name in sheet_names:
        sheet_df = pd.read_excel(metadata_xlsx, sheet_name=sheet_name, index_col=0)
        sheet_df = sheet_df.unstack().reset_index(name=sheet_name)  # unpivot (linearize) the table
        sheet_df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
        if plate_info_df.empty:
            plate_info_df = sheet_df
        else:
            plate_info_df = pd.merge(plate_info_df,
                                     sheet_df,
                                     how='left', on=['row_id', 'col_id'])
    antigen_df = antigen2D_to_df1D(xlsx_path=metadata_xlsx, sheet='antigen_array', data_col='antigen')
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
# python_df.dropna(inplace=True)
# Also update sera type to reflect their identity.
# posseralist=python_df_fix['Sera ID'].isin(['pos 1','pos 2','pos 3','pos 4'])
# python_df.loc[posseralist,'type'] = 'Diagnostic'


# In[171]:
python_df = python_df.append(df_scn)
python_df.replace([np.inf, -np.inf], np.nan, inplace=True)
python_df.dropna(subset=['OD'], inplace=True)


# %% Do the same for second path
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %% Read antigen and plate info
sheet_names = ['serum ID',
               'serum dilution',
               'serum cat',
               'serum type',
               'secondary ID',
               'secondary dilution']
plate_info_df2 = pd.DataFrame()
with pd.ExcelFile(metadata_path2) as metadata_xlsx:
    # get sheet names that are available in metadata
    sheet_names = list(set(metadata_xlsx.sheet_names).intersection(sheet_names))
    for sheet_name in sheet_names:
        sheet_df = pd.read_excel(metadata_xlsx, sheet_name=sheet_name, index_col=0)
        sheet_df = sheet_df.unstack().reset_index(name=sheet_name)  # unpivot (linearize) the table
        sheet_df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
        if plate_info_df2.empty:
            plate_info_df2 = sheet_df
        else:
            plate_info_df2 = pd.merge(plate_info_df2,
                                      sheet_df,
                                      how='left', on=['row_id', 'col_id'])
    antigen_df2 = antigen2D_to_df1D(xlsx_path=metadata_xlsx, sheet='antigen_array', data_col='antigen')
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
python_df2 = python_df2.append(df_scn)
python_df2.replace([np.inf, -np.inf], np.nan, inplace=True)
python_df2.dropna(subset=['OD'], inplace=True)

# %% Do the same for third path
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %% Read antigen and plate info
sheet_names = ['serum ID',
               'serum dilution',
               'serum cat',
               'serum type',
               'secondary ID',
               'secondary dilution']
plate_info_df3 = pd.DataFrame()
with pd.ExcelFile(metadata_path3) as metadata_xlsx:
    # get sheet names that are available in metadata
    sheet_names = list(set(metadata_xlsx.sheet_names).intersection(sheet_names))
    for sheet_name in sheet_names:
        sheet_df = pd.read_excel(metadata_xlsx, sheet_name=sheet_name, index_col=0)
        sheet_df = sheet_df.unstack().reset_index(name=sheet_name)  # unpivot (linearize) the table
        sheet_df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
        if plate_info_df3.empty:
            plate_info_df3 = sheet_df
        else:
            plate_info_df3 = pd.merge(plate_info_df3,
                                      sheet_df,
                                      how='left', on=['row_id', 'col_id'])
    antigen_df3 = antigen2D_to_df1D(xlsx_path=metadata_xlsx, sheet='antigen_array', data_col='antigen')
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
# python_df3.dropna(inplace=True)
# Also update sera type to reflect their identity.
# posseralist=python_df3_fix['Sera ID'].isin(['pos 1','pos 2','pos 3','pos 4'])
# python_df3.loc[posseralist,'type'] = 'Diagnostic'


# In[171]:
python_df3 = python_df3.append(df_scn)
python_df3.replace([np.inf, -np.inf], np.nan, inplace=True)
python_df3.dropna(subset=['OD'], inplace=True)

# %% Do the same for fourth path
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# %% Read antigen and plate info
sheet_names = ['serum ID',
               'serum dilution',
               'serum cat',
               'serum type',
               'secondary ID',
               'secondary dilution']
plate_info_df4 = pd.DataFrame()
with pd.ExcelFile(metadata_path4) as metadata_xlsx:
    # get sheet names that are available in metadata
    sheet_names = list(set(metadata_xlsx.sheet_names).intersection(sheet_names))
    for sheet_name in sheet_names:
        sheet_df = pd.read_excel(metadata_xlsx, sheet_name=sheet_name, index_col=0)
        sheet_df = sheet_df.unstack().reset_index(name=sheet_name)  # unpivot (linearize) the table
        sheet_df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
        if plate_info_df4.empty:
            plate_info_df4 = sheet_df
        else:
            plate_info_df4 = pd.merge(plate_info_df4,
                                      sheet_df,
                                      how='left', on=['row_id', 'col_id'])
    antigen_df4 = antigen2D_to_df1D(xlsx_path=metadata_xlsx, sheet='antigen_array', data_col='antigen')
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
python_df4 = python_df4.append(df_scn)
python_df4.replace([np.inf, -np.inf], np.nan, inplace=True)
python_df4.dropna(subset=['OD'], inplace=True)

#%% Remove failed wells for all dataframes
failed_wells = ['B6', 'B8', 'B10', 'B12','F4','H2','H4', 'H6', 'H9', 'H12']
for failed_well in failed_wells:
    python_df = python_df[python_df.well_id != failed_well]

wells2keep = ['B6', 'B8', 'B10', 'B12']
python_df22 = python_df2[(python_df2['well_id'].isin(wells2keep))]

failed_wells3 = ['A6']
for failedwell in failed_wells3:
    python_df3 = python_df3[python_df3.well_id == failedwell]

failed_wells4 = ['A6']
for failedwell in failed_wells4:
    python_df4 = python_df4[python_df4.well_id != failedwell]

#%% Concatenate dataframes
frames = [python_df, python_df22, python_df3, python_df4]
stitchedpython_df = pd.concat(frames)
# ## Fit curves to above plots

#%%
import scipy.optimize as optimization
import itertools

def fourPL(x, A, B, C, D):
    return ((A-D)/(1.0+((x/C)**(B))) + D)
   
def fit2df(df, model, pipeline):
    """fit model to x, y data in dataframe. 
    Return a dataframe with fit x, y for plotting
    """
    sera = df['serum ID'].unique()
    antigens = df['antigen'].unique()
    secondaries = df['secondary ID'].unique()

    keys = itertools.product(sera, antigens, secondaries)
    df_fit = pd.DataFrame(columns=df.columns)
    for serum, antigen, secondary in keys:
        print(serum, antigen)
        sec_dilu_df = df[(df['serum ID']==serum) &
                    (df['antigen']==antigen) &
                    (df['secondary ID'] == secondary) &
                    (df['pipeline']==pipeline)]
        sec_dilutions = sec_dilu_df['secondary dilution'].unique()
        for sec_dilution in sec_dilutions:
            sub_df = sec_dilu_df[(sec_dilu_df['secondary dilution'] == sec_dilution)].reset_index(drop=True)
            df_fit_temp = pd.DataFrame()
            guess = [0, 1, 5e-4, 1]
            xdata = sub_df['serum dilution'].to_numpy()
            ydata = sub_df['OD'].to_numpy()
            params, params_covariance = optimization.curve_fit(model, xdata, ydata, guess, bounds=(0, np.inf), maxfev=1e5)
            x_input = np.logspace(np.log10(np.min(xdata)), np.log10(np.max(xdata)), 50)
            y_fit = fourPL(x_input, *params)

            df_fit_temp['serum dilution'] = x_input
            df_fit_temp['OD'] = y_fit
            df_fit_temp['serum ID'] = ' '.join([serum, 'fit'])
            sub_df_expand = pd.concat(
                [sub_df.loc[[0], ['antigen',
                             'serum type', 'serum cat',
                             'secondary ID',
                             'secondary dilution',
                             'pipeline']]] * len(df_fit_temp.index), axis=0).reset_index(drop=True)
            df_fit_temp = pd.concat([df_fit_temp, sub_df_expand], axis=1)
            df_fit = df_fit.append(df_fit_temp)
    print('4PL fitting finished')
    return df_fit
python_df_fit = fit2df(stitchedpython_df, fourPL, 'python')
#%% plot the ODs and fits
# fig_path = os.path.join(data_folder1, 'pysero_plots3')
# os.makedirs(fig_path, exist_ok=True)
# # serum_type = 'Control'
# pipeline = 'python'
# #sera_list = natsorted(stitchedpython_df['serum ID'].unique())
# # 200612 JRB ADJUSTMENTS
#
# sera_list= ['Neg pool', 'Neg pool 1/200 + CR3022', 'Pos Pool']
# # 200612 JRB ADJUSTMENTS
# sera_fit_list = [' '.join([x, 'fit']) for x in sera_list]
# # sera_list = ['pos 1', 'pos 2', 'pos 3', 'pos 4', 'neg 1', 'neg 2', 'neg 3', 'neg 4']
# # sera_list = ['CR3022']
# # sera_list = ['Neg 1']
# # sera_list = ['neg ' + str(i) for i in range(11, 16)]
# # sera_list = ['positive ' + str(i) for i in range(1, 6)] + ['negative ' + str(i) for i in range(1, 7)]
# # markers = ['o'] * 5 + ['x'] * 6
# markers = 'o'
# # hue = 'secondary dilution'
# hue = 'serum ID'
# sec_dilutions = [2e-4]
# # sec_dilutions = [5e-5]
# # sec_dilutions = stitchedpython_df['secondary dilution'].unique()
# # style = 'secondary ID'
# style = 'serum type'
# antigens = natsorted(stitchedpython_df['antigen'].unique())
# serum_df = stitchedpython_df[(stitchedpython_df['pipeline']==pipeline) & stitchedpython_df['serum ID'].isin(sera_list)
#                      & stitchedpython_df['secondary dilution'].isin(sec_dilutions)]
# assert not serum_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
# for sec_id in serum_df['secondary ID'].unique():
#     sub_df = serum_df[(serum_df['secondary ID'] == sec_id)]
#     palette = sns.color_palette(n_colors=len(sub_df[hue].unique()))
#     print('plotting...')
#     g = sns.lmplot(x="serum dilution", y="OD", col_order=antigens,
#                     hue=hue, hue_order=sera_list, col="antigen", ci='sd', palette=palette, markers=markers,
#                      data=sub_df, col_wrap=5, fit_reg=False, x_estimator=np.mean)
#     sub_python_df_fit=python_df_fit[(python_df_fit['pipeline']==pipeline) &
#                                    python_df_fit['serum ID'].isin(sera_fit_list) &
#                                     (python_df_fit['secondary ID'] == sec_id) &
#                                     python_df_fit['secondary dilution'].isin(sec_dilutions)
#                                     ]
#     palette = sns.color_palette(n_colors=len(sub_python_df_fit[hue].unique()))
#     for antigen, ax in zip(antigens, g.axes.flat):
#         df_fit = sub_python_df_fit[(sub_python_df_fit['antigen'] == antigen)]
#         sns.lineplot(x="serum dilution", y="OD", hue=hue, hue_order=sera_fit_list, data=df_fit,
#                      style=style, palette=palette,
#                      ax=ax, legend=False)
#         ax.set(xscale="log")
#         # ax.set(ylim=[-0.05, 1.5])
#
#     plt.savefig(os.path.join(fig_path, '{}_{}_{}_fit.jpg'.format('THREE_sera', sec_id, sec_dilutions)),
#                              dpi=300, bbox_inches='tight')
#%% functions to compute ROC curves and AUC
# from sklearn.metrics import roc_curve, roc_auc_score
#
# def roc_from_df(df):
#     s = {}
#     y_test = df['serum type']
#     y_prob = df['OD']
#     s['False positive rate'], s['True positive rate'], _ = roc_curve(y_test, y_prob, pos_label='positive')
#     # s['AUC'] = roc_auc_score(y_test, y_prob)
#     s['AUC'] = [roc_auc_score(y_test, y_prob)] * len(s['False positive rate'])
#     return pd.Series(s)
#
# def get_roc_df(df):
#     """fit model to x, y data in dataframe.
#     Return a dataframe with fit x, y for plotting
#     """
#     df = df[df['serum type'].isin(['positive', 'negative'])]
#     roc_df = df[['antigen',
#                  'serum type',
#                  'secondary ID',
#                  'secondary dilution',
#                  'serum dilution',
#                  'OD',
#                  'pipeline']].drop_duplicates()
#     roc_df = roc_df.groupby(['antigen',
#                              'secondary ID',
#                              'secondary dilution',
#                              'serum dilution',
#                              'pipeline']).apply(roc_from_df)
#     # roc_df = roc_df.reset_index()
#     roc_df = roc_df.apply(pd.Series.explode).astype(float).reset_index()
#     return roc_df
# #%%
# serum_df = stitchedpython_df[~stitchedpython_df['serum ID'].str.contains('CR3022|Blank', regex=True)]
# print(serum_df['serum dilution'].unique())
# test_df = serum_df[(serum_df['serum dilution']==4.8828125e-06) & (serum_df['antigen']=='GFP foldon')]
# #%%
# roc_df = get_roc_df(stitchedpython_df)

#%% Plot one categorical scatter plot
# pipeline = 'python'
# dilution_list2 = ['0.02', '0.005', '0.00125']
# sec_dilutions = stitchedpython_df['secondary dilution'].unique()
# fig_path = os.path.join(data_folder1, 'pysero_plots3')
# os.makedirs(fig_path, exist_ok=True)
# antigens = natsorted(stitchedpython_df['antigen'].unique())
# serum_df = stitchedpython_df[(stitchedpython_df['pipeline']==pipeline) & stitchedpython_df['serum dilution'].isin(dilution_list2) & stitchedpython_df['secondary dilution'].isin(sec_dilutions)]
# assert not serum_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
#
# # Draw a categorical scatterplot to show each observation
# #palette=sns.palplot(sns.color_palette())
# palette1=sns.palplot(sns.color_palette())
# #palette2=sns.palplot(sns.color_palette("Reds",6))
# sns.catplot(x="serum dilution", y="OD", hue="serum ID", col_order=antigens, col="antigen", kind="point", join=False, palette=palette1,data=serum_df, col_wrap = 5, ci='sd')
# #sns.catplot(x="serum dilution", y="OD", hue="serum ID", col_order=antigens, col="antigen", kind="point", join=False, palette=palette2,data=serum_df2, col_wrap = 5, estimator=np.mean, ci='sd')
# plt.savefig(os.path.join(fig_path, '{}_{}_{}.jpg'.format('cat by dilution', 'check', 'check')),
#                               dpi=300, bbox_inches='tight')

#%% Make 2 plots, one for a sequential colour scheme of CR3022 and one for the other samples.
# Plot 1: CR3022 sequence
pipeline = 'python'
dilution_list2 = ['0.02', '0.005', '0.00125']
sera_1 = ['CR3022 1 ug/mL','CR3022 0.25 ug/mL','CR3022 0.0625 ug/mL','CR3022 0.01563 ug/mL','CR3022 0.00391 ug/mL','CR3022 0.00098 ug/mL','CR3022 0.00024 ug/mL','CR3022 0 ug/mL']
sec_dilutions = stitchedpython_df['secondary dilution'].unique()
fig_path = os.path.join(data_folder1, 'pysero_plots3')
os.makedirs(fig_path, exist_ok=True)
antigens = natsorted(stitchedpython_df['antigen'].unique())

serum_df1 = stitchedpython_df[(stitchedpython_df['pipeline']==pipeline) & stitchedpython_df['serum dilution'].isin(dilution_list2) & stitchedpython_df['secondary dilution'].isin(sec_dilutions) & stitchedpython_df['serum ID'].isin(sera_1)]
assert not serum_df1.empty, 'Plotting dataframe is empty. Please check the plotting keys'

# Draw a categorical scatterplot to show each observation
#palette=sns.palplot(sns.color_palette())
palette1=sns.color_palette("Greens", 8)
#palette2=sns.palplot(sns.color_palette("Reds",6))
g=sns.catplot(x="serum dilution", y="OD", hue="serum ID", hue_order=sera_1,col_order=antigens, col="antigen", kind="point", join=False, palette=palette1,data=serum_df1, col_wrap = 5, ci='sd')
g=g.set(ylim=[-0.05, 0.5])
#sns.catplot(x="serum dilution", y="OD", hue="serum ID", col_order=antigens, col="antigen", kind="point", join=False, palette=palette2,data=serum_df2, col_wrap = 5, estimator=np.mean, ci='sd')
plt.savefig(os.path.join(fig_path, '{}_{}_{}.jpg'.format('cat by dilution zoom','CR3022', 'all')),
                              dpi=300, bbox_inches='tight')
# Plot 2
pipeline = 'python'
dilution_list2 = ['0.02', '0.005', '0.00125']
sera_1 = ['Neg Pool', 'Pos Pool']
sec_dilutions = stitchedpython_df['secondary dilution'].unique()
fig_path = os.path.join(data_folder1, 'pysero_plots3')
os.makedirs(fig_path, exist_ok=True)
antigens = natsorted(stitchedpython_df['antigen'].unique())
serum_df1 = stitchedpython_df[(stitchedpython_df['pipeline']==pipeline) & stitchedpython_df['serum dilution'].isin(dilution_list2) & stitchedpython_df['secondary dilution'].isin(sec_dilutions) & stitchedpython_df['serum ID'].isin(sera_1)]
assert not serum_df1.empty, 'Plotting dataframe is empty. Please check the plotting keys'

# Draw a categorical scatterplot to show each observation
#palette=sns.palplot(sns.color_palette())
palette2=sns.color_palette()
#palette2=sns.palplot(sns.color_palette("Reds",6))
g=sns.catplot(x="serum dilution", y="OD", hue="serum ID", col_order=antigens, col="antigen", kind="point", join=False, palette=palette2,data=serum_df1, col_wrap = 5, ci='sd')
g.set(ylim=[-0.05, 0.5])
plt.savefig(os.path.join(fig_path, '{}_{}_{}.jpg'.format('cat by dilution zoom', 'pos pool', 'neg pool')),
                              dpi=300, bbox_inches='tight')
#%% Plot ROC curves
# hue = 'secondary dilution'
# hue = "serum dilution"
# # hue = 'serum ID'
# sec_dilutions = [2e-4]
# # sec_dilutions = [5e-5]
# # sec_dilutions = roc_df['secondary dilution'].unique()
# # style = 'secondary ID'
# style = 'serum type'
# pipeline='python'
# antigens = natsorted(roc_df['antigen'].unique())
# roc_py_df = roc_df[(roc_df['pipeline']==pipeline)
#                 & roc_df['secondary dilution'].isin(sec_dilutions)]
# sns.set_context("notebook")
# assert not roc_py_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
# for sec_id in roc_py_df['secondary ID'].unique():
#     palette = sns.color_palette(n_colors=len(sub_df[hue].unique()))
#     print('plotting...')
#
#     g = sns.FacetGrid(roc_py_df, hue=hue, col="antigen", col_wrap=5, aspect=1.05,
#                       hue_kws={'linestyle': ['-', '--', '-.', ':']})
#     g = (g.map(plt.plot, 'False positive rate', 'True positive rate').add_legend())
#
#     plt.savefig(os.path.join(fig_path, 'ROC_{}_{}.jpg'.format(sec_id, sec_dilutions)),
#                              dpi=300, bbox_inches='tight')


#%%
# #%% Generate plots from pysero
#
# # 4 positive sera vs 4 negative sera for control antigens
#
# fig_path = os.path.join(data_folder, 'pysero_fit_OD_8_sera_per_ag')
#
#
# sera_list = ['pos 1', 'pos 2', 'pos 3', 'pos 4', 'neg 1', 'neg 2', 'neg 3', 'neg 4'] # All +ve sera are affected by comets at some antigens, first 4 negative sera are not affected by comets.
# control_ag = ['xkappa-biotin','xIgG Fc', 'GFP foldon']
#
# # selectsera_ctlag = python_df_fix[(python_df_fix['Sera ID'].isin(seralist) & python_df_fix['antigen'].isin(controlag))]
# # selectsera_diagag = python_df_fix[(python_df_fix['Sera ID'].isin(seralist) & ~python_df_fix['antigen'].isin(controlag))]
#
# python_sera_df = python_df[python_df['serum ID'].isin(sera_list) & (python_df['pipeline']==pipeline)]
# y_list = ["OD", "intensity", "background"]
#
# for y in y_list:
#     g = sns.relplot(x="dilution", y=y,
#                         hue="serum ID", style = "type", col="antigen", ci=None,
#                          data=python_sera_df, col_wrap=3, estimator=np.nanmedian, kind='line')
#     g.set(xscale="log")
#     # g.set(xscale="log", ylim=[-0.05, 0.8], xlim=[1E-6, 2E-2])
#     plt.savefig(os.path.join(fig_path, '_'.join(['pysero', y + '.jpg'])), dpi=300, bbox_inches='tight')
#
#
#
# #%% Compare pysero and Scienion aross all serum and antigen combinations
# fig_path = os.path.join(data_folder, 'pysero_plots', 'pysero_metadata_vs_scienion')
# os.makedirs(fig_path, exist_ok=True)
# y_list = ["OD", "intensity", "background"]
# for serum in python_df['serum ID'].unique():
#     serum_df = python_df[python_df['serum ID']==serum]
#     for y in y_list:
#         g = sns.relplot(x="serum dilution", y=y,
#                     hue="pipeline", col="antigen", ci='sd', style='pipeline',
#                      data=serum_df, col_wrap=3, estimator=np.mean, kind='line')
#         g.set(xscale="log", ylim=[-0.05, 1.2])
#         plt.savefig(os.path.join(fig_path, '_'.join([serum, y + '.jpg'])), dpi=300, bbox_inches='tight')
#         plt.close()
#
#

