#!/usr/bin/env python
# coding: utf-8

#%% Setup


from pprint import pprint
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib
import skimage.io as io
import seaborn as sns; sns.set_context("talk")
font = {'size'   : 16, 'weight' : 'normal', 'family' : 'arial'}
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
    df = df.unstack().reset_index(name=data_col) # unpivot (linearize) the table
    df.rename(columns={'level_1': 'antigen_row', 'level_0': 'antigen_col'}, inplace=True)
    df[['antigen_row', 'antigen_col']] = df[['antigen_row', 'antigen_col']].applymap(int)
    df = df[['antigen_row', 'antigen_col', data_col]]
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
    df = df.unstack().reset_index(name=data_col) # unpivot (linearize) the table
    df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
    df['well_id'] = df.row_id + df.col_id.map(str)
    df = df[['well_id', data_col]]
    return df
#%% Set paths
data_folder=r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-04-08-14-48-49-COVID_concentrationtest_April8_images'
metadata_path=os.path.join(data_folder, 'pysero_output_data_metadata.xlsx')
OD_path=os.path.join(data_folder, 'median_ODs.xlsx')
int_path=os.path.join(data_folder, 'median_intensities.xlsx')
bg_path=os.path.join(data_folder, 'median_backgrounds.xlsx')
# antigens_path=os.path.join(data_folder,'python_antigens.xlsx')
scienion_path=os.path.join(data_folder, 'Scienion_reader_output', '2020-04-08-14-57-09-COVID_concentrationtest_April8_analysis2.xlsx')
#%% Read antigen and plate info
sheet_names = ['serum ID', 'serum dilution', 'serum type']
plate_info_df = pd.DataFrame()
for sheet_name in sheet_names:
    
    sheet_df = pd.read_excel(metadata_path, sheet_name=sheet_name, index_col=0)
    sheet_df = sheet_df.unstack().reset_index(name=sheet_name) # unpivot (linearize) the table
    sheet_df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
    if sheet_name == 'serum ID':
        plate_info_df = sheet_df
    else:
        plate_info_df = pd.merge(plate_info_df,
                                 sheet_df,
                                 how='left', on=['row_id', 'col_id'])       
plate_info_df['well_id'] = plate_info_df.row_id + plate_info_df.col_id.map(str)
plate_info_df['serum dilution'] = 1 / plate_info_df['serum dilution']
plate_info_df = plate_info_df[['serum ID', 'serum dilution', 'serum type', 'well_id']]
plate_info_df.dropna(inplace=True)
# # serum_day_df=plate_info_df['Sera ID'].str.split('-Day', expand=True) #If sera names are not organized by ID+collection day, comment out
# #serum_day_df.fillna(0, inplace=True)
# #serum_day_df.columns=['serum', 'day']
# #serum_day_df['day'] = serum_day_df['day'].astype(int)
# #plate_info_df = pd.concat([plate_info_df, serum_day_df], axis=1)
# #plate_info_df.drop(['Sera ID'], axis=1, inplace=True)
# plate_info_df['type'] = plate_info_df['Sera ID'].str.replace(' \d+', '')
# # print(plate_info_df[plate_info_df['Sera ID']=='pos 1']['dilution'].unique())
#%% Read antigen information.
antigen_df = antigen2D_to_df1D(xlsx_path=metadata_path, sheet='antigen_array', data_col='antigen')
#%% Read analysis output from Scienion
# Read all wells into dictionary. 
scienion_df = pd.DataFrame()
for well_id in plate_info_df['well_id']:
    OD_1_antiten_df = pd.read_excel(scienion_path, sheet_name=well_id)
    OD_1_antiten_df['well_id'] = well_id
    scienion_df = scienion_df.append(OD_1_antiten_df, ignore_index=True)
#%%   parse spot ids
spot_id_df=scienion_df['ID'].str.extract(r'spot-(\d)-(\d)')
spot_id_df = spot_id_df.astype(int) - 1 # index starting from 0
spot_id_df.rename(columns={0: 'antigen_row', 1: 'antigen_col'}, inplace=True)
scienion_df = pd.concat([spot_id_df, scienion_df], axis=1)
scienion_df.drop('ID', axis=1, inplace=True)

#%% invert the intensity and compute ODs, check A2
df_scn = scienion_df.loc[:, ['antigen_row', 'antigen_col', 'well_id']]
df_scn['intensity'] = 1 - scienion_df['Median'] / 255
df_scn['background'] = 1 - scienion_df['Background Median'] / 255
df_scn['OD'] = np.log10(df_scn['background'] / df_scn['intensity'])
#%% Join Scienion data with plateInfo


df_scn = pd.merge(df_scn,
                 antigen_df,
                 how='left', on=['antigen_row', 'antigen_col'])
df_scn = pd.merge(df_scn,
                 plate_info_df,
                 how='right', on=['well_id'])


#%% Add a name of the pipeline to the dataframe.

df_scn['pipeline'] = 'scienion'


#%% Read optical density from pysero
OD_df = pd.DataFrame()
int_df = pd.DataFrame()
bg_df = pd.DataFrame()
for _, row in antigen_df.iterrows():
    sheet_name = 'od_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
    OD_1_antiten_df = well2D_to_df1D(xlsx_path=OD_path, sheet=sheet_name, data_col='OD')
    OD_1_antiten_df['antigen_row'] = row['antigen_row']
    OD_1_antiten_df['antigen_col'] = row['antigen_col']
    OD_1_antiten_df['antigen'] = row['antigen']
    OD_df = OD_df.append(OD_1_antiten_df, ignore_index=True)

for _, row in antigen_df.iterrows():
    sheet_name = 'int_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
    int_1_antiten_df = well2D_to_df1D(xlsx_path=int_path, sheet=sheet_name, data_col='intensity')
    int_1_antiten_df['antigen_row'] = row['antigen_row']
    int_1_antiten_df['antigen_col'] = row['antigen_col']
    int_1_antiten_df['antigen'] = row['antigen']
    int_df = int_df.append(int_1_antiten_df, ignore_index=True)

for _, row in antigen_df.iterrows():
    sheet_name = 'bg_{}_{}_{}'.format(row['antigen_row'], row['antigen_col'], row['antigen'])
    bg_1_antiten_df = well2D_to_df1D(xlsx_path=bg_path, sheet=sheet_name, data_col='bg')
    bg_1_antiten_df['antigen_row'] = row['antigen_row']
    bg_1_antiten_df['antigen_col'] = row['antigen_col']
    bg_1_antiten_df['antigen'] = row['antigen']
    bg_df = bg_df.append(bg_1_antiten_df, ignore_index=True)

#%% merge OD with antigen and plate info.


# Use of filter avoids merge of duplicate columns when the cell is run multiple times.
OD_df = OD_df.filter(items=['antigen_row', 'antigen_col','OD','well_id'],axis=1)
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

#%% Generate plots from pysero

# 4 positive sera vs 4 negative sera for control antigens


fig_path = os.path.join(data_folder, 'pysero_fit_OD_8_sera_per_ag')
os.makedirs(fig_path, exist_ok=True)
pipeline = 'python'
sera_list = ['pos 1', 'pos 2', 'pos 3', 'pos 4', 'neg 1', 'neg 2', 'neg 3', 'neg 4'] # All +ve sera are affected by comets at some antigens, first 4 negative sera are not affected by comets.
control_ag = ['xkappa-biotin','xIgG Fc', 'GFP foldon']

# selectsera_ctlag = python_df_fix[(python_df_fix['Sera ID'].isin(seralist) & python_df_fix['antigen'].isin(controlag))]
# selectsera_diagag = python_df_fix[(python_df_fix['Sera ID'].isin(seralist) & ~python_df_fix['antigen'].isin(controlag))]

python_sera_df = python_df[python_df['serum ID'].isin(sera_list) & (python_df['pipeline']==pipeline)]
y_list = ["OD", "intensity", "background"]

for y in y_list:
    g = sns.relplot(x="dilution", y=y,
                        hue="serum ID", style = "type", col="antigen", ci=None, 
                         data=python_sera_df, col_wrap=3, estimator=np.nanmedian, kind='line')
    g.set(xscale="log")
    # g.set(xscale="log", ylim=[-0.05, 0.8], xlim=[1E-6, 2E-2])
    plt.savefig(os.path.join(fig_path, '_'.join(['pysero', y + '.jpg'])), dpi=300, bbox_inches='tight')



#%% Compare pysero and Scienion aross all serum and antigen combinations


fig_path = os.path.join(data_folder, 'pysero_plots', 'pysero_metadata_vs_scienion')
os.makedirs(fig_path, exist_ok=True)
y_list = ["OD", "intensity", "background"]
for serum in python_df['serum ID'].unique():
    serum_df = python_df[python_df['serum ID']==serum]
    for y in y_list:
        g = sns.relplot(x="serum dilution", y=y,
                    hue="pipeline", col="antigen", ci='sd', style='pipeline',
                     data=serum_df, col_wrap=3, estimator=np.mean, kind='line')
        g.set(xscale="log", ylim=[-0.05, 1.2])
        plt.savefig(os.path.join(fig_path, '_'.join([serum, y + '.jpg'])), dpi=300, bbox_inches='tight')
        plt.close()


# ## Fit curves to above plots

# In[154]:


# Compare sera for each antigen
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
    keys = itertools.product(sera, antigens)
    df_fit = pd.DataFrame(columns=df.columns)
    for serum, antigen in keys:
        print(serum, antigen)
        sub_df = df[(df['serum ID']==serum) & (df['antigen']==antigen) & (df['pipeline']==pipeline)]
        df_fit_temp = pd.DataFrame(columns=df.columns)
        guess = [0, 1, 5e-4, 1]
        xdata = sub_df['dilution'].to_numpy()
        ydata = sub_df['OD'].to_numpy()
        params, params_covariance = optimization.curve_fit(model, xdata, ydata, guess, bounds=(0, np.inf), maxfev=1e5)
        x_input = np.logspace(np.log10(np.min(xdata)), np.log10(np.max(xdata)), 50)
        y_fit = fourPL(x_input, *params)

        df_fit_temp['dilution'] = x_input
        df_fit_temp['OD'] = y_fit
        df_fit_temp['serum ID'] = ' '.join([serum, 'fit'])
        df_fit_temp['antigen'] = antigen
        df_fit_temp['type'] = sub_df['type'].unique()[0]
        df_fit_temp['pipeline'] = pipeline
        df_fit = df_fit.append(df_fit_temp)  
    return df_fit
python_df_fit = fit2df(python_df, fourPL, 'python')


# In[157]:


# serum_type = 'Control'
sera_list = ['pos 1', 'pos 2', 'pos 3', 'pos 4', 'neg 1', 'neg 2', 'neg 3', 'neg 4'] 
# sera_list = ['neg ' + str(i) for i in range(11, 16)]
sub_df = python_df[(python_df['pipeline']==pipeline) &
                  python_df['serum ID'].isin(sera_list) ]

g = sns.lmplot(x="dilution", y="OD",
                hue='serum ID', col="antigen", ci='sd',
                 data=sub_df, col_wrap=3, fit_reg=False, x_estimator=np.mean)

sera_list = [' '.join([x, 'fit']) for x in sera_list]
sub_python_df_fit=python_df_fit[(python_df_fit['pipeline']==pipeline) & 
                               python_df_fit['serum ID'].isin(sera_list) ]

for antigen, ax in zip(sub_df['antigen'].unique(), g.axes.flat):
    df_fit = sub_python_df_fit[sub_python_df_fit['antigen']==antigen]
    sns.lineplot(x="dilution", y="OD", hue='serum ID', data=df_fit, style='type', ax=ax, legend=False)
    ax.set(xscale="log", ylim=[-0.05, 1.5])
    
# plt.savefig(os.path.join(fig_path, 'pyseroOD_neg_sera_11_16_fit.jpg'), dpi=300, bbox_inches='tight')
plt.savefig(os.path.join(fig_path, 'pyseroOD_4_pos_4_neg_sera_fit.jpg'), dpi=300, bbox_inches='tight')





