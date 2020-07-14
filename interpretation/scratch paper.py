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
import scipy.optimize as optimization
import itertools

import math
d = {'OD': [0.9, 0.5,0.2,0.1,0.99,0.01], 'antigen': ['RBD','Spike','RBD','Spike','RBD','Spike'],'serum ID1':['test1','test1','test1','test1','Pos Pool','Pos Pool'],'sc':['test1','test1','test1','test1','test1','test1'],'value on curve':[0.3,0.4,0.1,0.4,'',0.1]}
first_df = pd.DataFrame(data=d)
print(first_df)

rangeforag_df = first_df[(first_df['serum ID1'] == first_df['sc'])]
minmaxValue_df = pd.DataFrame()
for antigen in rangeforag_df['antigen'].unique():
    antigen_df = rangeforag_df[(rangeforag_df['antigen'] ==antigen)]
    maxValue = antigen_df['OD'].max()
    minValue = antigen_df['OD'].min()
    print(maxValue)
    temp_df = pd.DataFrame()
    temp_df['maxValue'] = [maxValue]
    temp_df['minValue'] = [minValue]
    temp_df['antigen'] = [antigen]

    minmaxValue_df = minmaxValue_df.append(temp_df)
print(minmaxValue_df)
first_df = first_df.merge(minmaxValue_df, left_on='antigen', right_on='antigen', how='outer')
print(first_df)

#if OD is greater than max value, change value on curve to 1
# Do the converse for min value, change to 0.000000001
for idx, row in first_df.iterrows():
    if row['OD'] > row['maxValue']:
        # first_df.loc((row['value on curve']),1)
        first_df.loc[idx,'value on curve'] = 1
    elif row['OD'] < row['minValue']:

        first_df.loc[idx, 'value on curve'] = 0.000000001


print(first_df)


