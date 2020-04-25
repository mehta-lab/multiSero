#!/usr/bin/env python
# coding: utf-8

# In[1]:


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

def xlsx2D_to_df1D(xlsx_path, sheet, data_col):
    df = pd.read_excel(xlsx_path, sheet_name=sheet, index_col=0)
    df = df.unstack().reset_index(name=data_col) # unpivot (linearize) the table
    df.rename(columns={'level_1': 'antigen_row', 'level_0': 'antigen_col'}, inplace=True)
    df = df[['antigen_row', 'antigen_col', data_col]]
    return df
    


# # Setup

# *  xkappa biotin spots: these spots are coated with human kappa light chain conjugated to biotin, so that when the streptavidin-HRP is added there is strong signal independent of the serum content.  **We should not expect the fiducial signal to weaken with serial dilution.**
# *  Row-5 (row index 4) in multiple wells shows excessive comets and are not interpretable. The row corresponds to concentration range of (SARS-CoV2-RBD).  
# *  Looking at debug plots, following spots are not reliable:
#  * A1: row-4, (5,1)
#  * B1: row-4, (5,1)
#  * C3: whole well
#  * D3: whole well
#  * E3: whole well
#  * E11: whole well
# 

# ## Set paths

# In[108]:


data_folder=r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-04-08-14-48-49-COVID_concentrationtest_April8_images'
plateinfo_path=os.path.join(data_folder,'COVIDplateinfo0408_v2.xlsx')
antigenOD_path=os.path.join(data_folder,'python_fit_ODs.xlsx') 
python_int_path=os.path.join(data_folder,'python_fit_intensities.xlsx') 
python_bg_path=os.path.join(data_folder,'python_fit_backgrounds.xlsx')
antigens_path=os.path.join(data_folder,'python_antigens.xlsx')
scienionpath=os.path.join(data_folder,'2020-04-08-14-57-09-COVID_concentrationtest_April8_analysis2.xlsx')


# # Read antigen and plate info

# In[161]:



sheet_names = ['serum ID', 'dilution', 'type']
plate_info_df = pd.DataFrame()
for sheet_name in sheet_names:
    
    sheet_df = pd.read_excel(plateinfo_path, sheet_name=sheet_name, index_col=0)
    sheet_df = sheet_df.unstack().reset_index(name=sheet_name) # unpivot (linearize) the table
    sheet_df.rename(columns={'level_1': 'row_id', 'level_0': 'col_id'}, inplace=True)
    if sheet_name == 'serum ID':
        plate_info_df = sheet_df
    else:
        plate_info_df = pd.merge(plate_info_df,
                                 sheet_df,
                                 how='left', on=['row_id', 'col_id'])       
plate_info_df['well_id'] = plate_info_df.row_id + plate_info_df.col_id.map(str)
plate_info_df['dilution'] = 1 / plate_info_df['dilution']  
plate_info_df = plate_info_df[['serum ID', 'dilution', 'type', 'well_id']]
plate_info_df.dropna(inplace=True)
# plate_info_df.rename(columns={'Well': 'well_id', 'Type':'type', 'Dilution':'dilution'}, inplace=True)

# # serum_day_df=plate_info_df['Sera ID'].str.split('-Day', expand=True) #If sera names are not organized by ID+collection day, comment out
# #serum_day_df.fillna(0, inplace=True)
# #serum_day_df.columns=['serum', 'day']
# #serum_day_df['day'] = serum_day_df['day'].astype(int)
# #plate_info_df = pd.concat([plate_info_df, serum_day_df], axis=1)
# #plate_info_df.drop(['Sera ID'], axis=1, inplace=True)
# plate_info_df['type'] = plate_info_df['Sera ID'].str.replace(' \d+', '')
# # print(plate_info_df[plate_info_df['Sera ID']=='pos 1']['dilution'].unique())



# In[117]:


antigen_df = xlsx2D_to_df1D(xlsx_path=antigens_path, sheet='antigens', data_col='antigen')   



# # Read analysis output from Scienion 

# In[123]:


# Read all wells into dictionary. 
scienion_df = pd.DataFrame()
for well_id in plate_info_df['well_id']:
    OD_1_well_df = pd.read_excel(scienionpath, sheet_name=well_id)
    OD_1_well_df['well_id'] = well_id
    scienion_df = scienion_df.append(OD_1_well_df, ignore_index=True)



# In[124]:


# parse spot ids
spot_id_df=scienion_df['ID'].str.extract(r'spot-(\d)-(\d)')
spot_id_df = spot_id_df.astype(int) - 1 # index starting from 0
spot_id_df.rename(columns={0: 'antigen_row', 1: 'antigen_col'}, inplace=True)

scienion_df = pd.concat([spot_id_df, scienion_df], axis=1)
scienion_df.drop('ID', axis=1, inplace=True)


# In[162]:


# invert the intensity and compute ODs, check A2
df_scn = scienion_df.loc[:, ['antigen_row', 'antigen_col', 'well_id']]
df_scn['intensity'] = 1 - scienion_df['Median'] / 255
df_scn['background'] = 1 - scienion_df['Background Median'] / 255
df_scn['OD'] = np.log10(df_scn['background'] / df_scn['intensity'])
df_scn[df_scn['well_id']=='A1']


# ### Join Scienion data with plateInfo

# In[163]:


df_scn = pd.merge(df_scn,
                 antigen_df,
                 how='left', on=['antigen_row', 'antigen_col'])
df_scn = pd.merge(df_scn,
                 plate_info_df,
                 how='right', on=['well_id'])


# In[165]:


df_scn['pipeline'] = 'scienion'
df_scn


# # Read output from pysero

# ## Read optical density

# In[128]:


# Read all wells into dictionary and into a 4D numpy array.
OD_df = pd.DataFrame()
for well_id in plate_info_df['well_id']:
    OD_1_well_df = xlsx2D_to_df1D(xlsx_path=antigenOD_path, sheet=well_id, data_col='OD')   
    OD_1_well_df['well_id'] = well_id
    OD_df = OD_df.append(OD_1_well_df, ignore_index=True)
OD_df


# ### merge OD with antigen and plate info.

# In[167]:


# Use of filter avoids merge of duplicate columns when the cell is run multiple times.
OD_df = OD_df.filter(items=['antigen_row', 'antigen_col','OD','well_id'],axis=1)
OD_df = pd.merge(OD_df,
                 antigen_df,
                 how='left', on=['antigen_row', 'antigen_col'])
OD_df = pd.merge(OD_df,
                 plate_info_df,
                 how='right', on=['well_id'])


# ## median intensities and background from pysero

# In[168]:


# Read all intensity wells into dictionary and into a 4D numpy array.
python_int_df = pd.DataFrame()
for well_id in plate_info_df['well_id']:
    OD_1_well_df = xlsx2D_to_df1D(xlsx_path=python_int_path, sheet=well_id, data_col='intensity')   
    OD_1_well_df['well_id'] = well_id
    python_int_df = python_int_df.append(OD_1_well_df, ignore_index=True)
python_int_df


# In[169]:


# Read all background wells into dictionary and into a 4D numpy array.
python_bg_df = pd.DataFrame()
for well_id in plate_info_df['well_id']:
    OD_1_well_df = xlsx2D_to_df1D(xlsx_path=python_bg_path, sheet=well_id, data_col='background')   
    OD_1_well_df['well_id'] = well_id
    python_bg_df = python_bg_df.append(OD_1_well_df, ignore_index=True)
python_bg_df


# ### merge intensity and background with OD, plateinfo, antigen info

# In[170]:


# OD_df = OD_df.filter(items=['well_id','antigen_row', 'antigen_col','antigen','Sera ID', 'type', 'dilution', 'OD'],axis=1)
python_df = pd.merge(OD_df,
                 python_int_df,
                 how='left', on=['antigen_row', 'antigen_col', 'well_id'])
python_df = pd.merge(python_df,
                 python_bg_df,
                 how='left', on=['antigen_row', 'antigen_col', 'well_id'])

python_df['pipeline'] = 'python'
# python_df.dropna(inplace=True)
# Also update sera type to reflect their identity.
# posseralist=python_df_fix['Sera ID'].isin(['pos 1','pos 2','pos 3','pos 4'])
# python_df.loc[posseralist,'type'] = 'Diagnostic'


# In[171]:


python_df = python_df.append(df_scn)


# In[172]:


python_df


# # Generate plots from pysero

# ## 4 positive sera vs 4 negative sera for control antigens

# In[148]:


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


# ## Compare pysero and Scienion aross all serum and antigen combinations

# In[180]:


fig_path = os.path.join(data_folder, 'python_median_vs_scienion')
os.makedirs(fig_path, exist_ok=True)
y_list = ["OD", "intensity", "background"]
for serum in python_df['serum ID'].unique():
    serum_df = python_df[python_df['serum ID']==serum]
    for y in y_list:
        g = sns.relplot(x="dilution", y=y,
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


# In[142]:


python_df_fit


# In[ ]:




