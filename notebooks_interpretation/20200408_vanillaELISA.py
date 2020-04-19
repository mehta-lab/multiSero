#!/usr/bin/env python
# coding: utf-8

# In[]:

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
    df.rename(columns={'level_0': 'well_col', 'level_1': 'well_row'}, inplace=True)
    df['well_col'] = df['well_col'].apply(str)
    df['well_id'] =  df['well_row']+df['well_col']
    df = df[['well_id', data_col]]
    return df

# Setup

# * RBD is coated at 2ug/ml for each well
# * The positive samples are serial diluted down each column
# * The secondary antibody are serial diluted across rows
# 

# ## Set paths

# In[15]:


data_folder=r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_leica/20200408_VanillaELISA_Leica2.5x_2'
plateinfo_path=os.path.join(data_folder,'VanillaElisaPlateInfo.csv')
#antigenOD_path=os.path.join(data_folder,'python_median_ODs.xlsx') 
python_int_path=os.path.join(data_folder,'intensities.xlsx')
#python_bg_path=os.path.join(data_folder,'python_median_backgrounds.xlsx')
#antigens_path=os.path.join(data_folder,'python_antigens.xlsx')
#scienionpath=os.path.join(data_folder,'2020-04-08-14-57-09-COVID_concentrationtest_April8_analysis2.xlsx')


# # Read antigen and plate info

# In[18]:


df = xlsx2D_to_df1D(python_int_path, 'od', 'od')

# plate_info_df=pd.read_csv(plateinfo_path)
# plate_info_df.rename(columns={'Well': 'well_id','Antibody ID':'antibody_id' ,'Type':'type', 'Antibody Dilution':'antibody_dilution','Sample ID ':'sample_id','Sample Concentration':'sample_dilution', 'Intensity':'intensity'}, inplace=True)

#serum_day_df=plate_info_df['Sera ID'].str.split('-Day', expand=True) #If sera names are not organized by ID+collection day, comment out
#serum_day_df.fillna(0, inplace=True)
#serum_day_df.columns=['serum', 'day']
#serum_day_df['day'] = serum_day_df['day'].astype(int)
#plate_info_df = pd.concat([plate_info_df, serum_day_df], axis=1)
#plate_info_df.drop(['Sera ID'], axis=1, inplace=True)

#print(plate_info_df[plate_info_df['sample_id']=='Sample 1']['antibody_dilution'].unique())
#plate_info_df


# In[7]:


#This plate has 1 antigen, RBD coated at 2 ug/mL
#antigen_df = xlsx2D_to_df1D(xlsx_path=antigens_path, sheet='antigens', data_col='antigen')   
#antigen_df


# # Generate plots

# ## Plot Antibody dilution series per sample

# In[41]:


fig_path = os.path.join(data_folder, 'antibody_dilutions_per_sample')
os.makedirs(fig_path, exist_ok=True)

samplelist = ['Sample 1','Sample 2']
controlsamples = ['Blank','Negative']

selectsamples_exps = plate_info_df[(plate_info_df['sample_id'].isin(samplelist))]
selectsamples_ctls = plate_info_df[(plate_info_df['sample_id'].isin(controlsamples))]

g = sns.relplot(x="sample_dilution", y="intensity", hue = "antibody_dilution", col="sample_id", ci=None, data=selectsamples_exps, col_wrap=2, estimator=np.nanmedian, kind='line')
#g.set(xscale="log", ylim=[6E6, 65E6], xlim=[1E-6, 0.01])
g.set(xscale="log", xlim=[1E-4, 0.004])


#g = sns.relplot(x="sample_dilution", y="intensity", hue = "antibody_dilution", col="sample_id", ci=None, data=selectsamples_ctls, col_wrap=2, estimator=np.nanmedian, kind='line')
#g.set(xscale="log", xlim=[0.01, 1], ylim=[13000,63000])
plt.savefig(os.path.join(fig_path, '_'.join(['antibody_dilutions_per_sample' + '.jpg'])), dpi=300, bbox_inches='tight')


# ## Plot sample dilution series per secondary antibody

# In[45]:


fig_path = os.path.join(data_folder, 'sample_dilutions_per_antibody')
os.makedirs(fig_path, exist_ok=True)

g = sns.relplot(x="antibody_dilution", y="intensity", hue = "sample_dilution", col="antibody_id", ci=None, data=selectsamples_exps, col_wrap=2, estimator=np.nanmedian, kind='line')
#g.set(xscale="log", ylim=[6E6, 65E6], xlim=[1E-6, 0.01])
g.set(xscale="log", xlim=[1E-6, 0.004])

g = sns.relplot(x="antibody_dilution", y="intensity", hue = "sample_dilution", col="antibody_id", ci=None, data=selectsamples_ctls, col_wrap=2, estimator=np.nanmedian, kind='line')
g.set(xscale="log", xlim=[1E-6, 0.004])

plt.savefig(os.path.join(fig_path, '_'.join(['sample_dilutons_per_antibody' + '.jpg'])), dpi=300, bbox_inches='tight')


# ## Fit curves to above plots

# In[50]:


g = sns.lmplot(x="antibody_dilution", y="intensity", hue = "sample_dilution", col="antibody_id", ci=None, data=selectsamples_exps, col_wrap=2, logistic=True,truncate=True)
g.set(xscale="log", ylim=[-0.01,0.8], xlim=[1E-6, 1E-2]);
plt.savefig(os.path.join(fig_path, '_'.join(['sampledilutionsperabFit' + '.jpg'])), dpi=300, bbox_inches='tight')


# In[162]:


g = sns.lmplot(x="dilution", y="OD", hue="Sera ID", col="antigen", col_wrap=3, data=selectsera_ctlag, ci=None, 
               logistic=True, truncate=True)
g.set(xscale="log", ylim=[-0.01,0.8], xlim=[1E-6, 1E-2]);


# In[166]:


g = sns.lmplot(x="dilution", y="OD", hue="Sera ID", col="antigen", col_wrap=3, data=selectsera_diagag, ci=None, 
                 robust=True, truncate = True)
g.set(xscale="log", ylim=[-0.01,0.8], xlim=[1E-6, 1E-2]);
plt.savefig(os.path.join(fig_path, '_'.join(['pyseroOD_DiagAgFit' + '.jpg'])), dpi=300, bbox_inches='tight')


# In[ ]:




