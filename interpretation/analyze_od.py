#!/usr/bin/env python
# coding: utf-8

# %% Setup


from pprint import pprint
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib
from natsort import natsorted
import seaborn as sns
import scipy.optimize as optimization
import itertools
from sklearn.metrics import roc_curve, roc_auc_score


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

def read_plate_info(metadata_xlsx):
    print('Reading the plate info...')

    sheet_names = ['serum ID',
                   'serum dilution',
                   'serum type',
                   'secondary ID',
                   'secondary dilution']
    plate_info_df = pd.DataFrame()
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
    plate_info_df['well_id'] = plate_info_df.row_id + plate_info_df.col_id.map(str)
    sheet_names.append('well_id')
    # convert to number and non-numeric to NaN
    plate_info_df['serum dilution'] = \
        plate_info_df['serum dilution'].apply(pd.to_numeric, errors='coerce')
    plate_info_df.dropna(inplace=True)
    if np.all(plate_info_df['serum dilution'] >= 1):
        # convert dilution to concentration
        plate_info_df['serum dilution'] = 1 / plate_info_df['serum dilution']
    plate_info_df.drop(['row_id', 'col_id'], axis=1, inplace=True)
    return plate_info_df

def read_antigen_info(metadata_path):
    print('Reading antigen information...')
    antigen_df = antigen2D_to_df1D(xlsx_path=metadata_path, sheet='antigen_array', data_col='antigen')
    antigen_type_df = antigen2D_to_df1D(xlsx_path=metadata_path, sheet='antigen_type', data_col='antigen type')
    antigen_df = pd.merge(antigen_df, antigen_type_df, how='left', on=['antigen_row', 'antigen_col'])
    return antigen_df

def read_pysero_output(file_path, antigen_df, file_type='od'):
    print('Reading {}...'.format(file_type))
    data_col = {'od': 'OD', 'int': 'intensity', 'bg': 'background'}
    data_df = pd.DataFrame()

    with pd.ExcelFile(file_path) as file:
        for _, row in antigen_df.iterrows():
            sheet_name = '{}_{}_{}_{}'.format(file_type, row['antigen_row'], row['antigen_col'], row['antigen'])
            data_1_antiten_df = well2D_to_df1D(xlsx_path=file, sheet=sheet_name, data_col=data_col[file_type])
            data_1_antiten_df['antigen_row'] = row['antigen_row']
            data_1_antiten_df['antigen_col'] = row['antigen_col']
            data_1_antiten_df['antigen'] = row['antigen']
            data_df = data_df.append(data_1_antiten_df, ignore_index=True)
    return data_df

def read_scn_output(file_path, plate_info_df):
    # Read analysis output from Scienion
    scienion_df = pd.DataFrame()
    with pd.ExcelFile(file_path) as scienion_xlsx:
        for well_id in plate_info_df['well_id']:
            OD_1_antiten_df = pd.read_excel(scienion_xlsx, sheet_name=well_id)
            OD_1_antiten_df['well_id'] = well_id
            scienion_df = scienion_df.append(OD_1_antiten_df, ignore_index=True)
    # parse spot ids
    spot_id_df = scienion_df['ID'].str.extract(r'spot-(\d)-(\d)')
    spot_id_df = spot_id_df.astype(int) - 1  # index starting from 0
    spot_id_df.rename(columns={0: 'antigen_row', 1: 'antigen_col'}, inplace=True)
    scienion_df = pd.concat([spot_id_df, scienion_df], axis=1)
    scienion_df.drop('ID', axis=1, inplace=True)
    # invert the intensity and compute ODs
    df_scn = scienion_df.loc[:, ['antigen_row', 'antigen_col', 'well_id']]
    df_scn['intensity'] = 1 - scienion_df['Median'] / 255
    df_scn['background'] = 1 - scienion_df['Background Median'] / 255
    df_scn['OD'] = np.log10(df_scn['background'] / df_scn['intensity'])
    return df_scn

def slice_df(df, slice_action, column, keys):
    if slice_action is None:
        return df
    elif slice_action == 'keep':
        df = df[df[column].isin(keys)]
    elif slice_action == 'drop':
        df = df[~df[column].isin(keys)]
    else:
        raise ValueError('slice action has to be "keep" or "drop", not "{}"'.format(slice_action))
    return df

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
                             'serum type',
                             'secondary ID',
                             'secondary dilution',
                             'pipeline']]] * len(df_fit_temp.index), axis=0).reset_index(drop=True)
            df_fit_temp = pd.concat([df_fit_temp, sub_df_expand], axis=1)
            df_fit = df_fit.append(df_fit_temp)
    print('4PL fitting finished')
    return df_fit

def roc_from_df(df):
    s = {}
    y_test = df['serum type']
    y_prob = df['OD']
    s['False positive rate'], s['True positive rate'], s['threshold'] = roc_curve(y_test, y_prob, pos_label='positive')
    try:
        s['AUC'] = [roc_auc_score(y_test, y_prob)] * len(s['False positive rate'])
    except ValueError as err:
        print('serum dilution {} only has {} serum type. {}'.
              format(df['serum dilution'].unique()[0], y_test.unique()[0], err))
        s['AUC'] = [np.nan] * len(s['False positive rate'])
    return pd.Series(s)

def get_roc_df(df):
    """fit model to x, y data in dataframe.
    Return a dataframe with fit x, y for plotting
    """
    df = df[df['serum type'].isin(['positive', 'negative'])]
    roc_df = df[['antigen',
                 'serum type',
                 'secondary ID',
                 'secondary dilution',
                 'serum dilution',
                 'OD',
                 'pipeline']]
    roc_df = roc_df.groupby(['antigen',
                             'secondary ID',
                             'secondary dilution',
                             'serum dilution',
                             'pipeline']).apply(roc_from_df)
    # roc_df = roc_df.reset_index()
    roc_df = roc_df.apply(pd.Series.explode).astype(float).reset_index()
    roc_df.dropna(inplace=True)
    return roc_df

def roc_plot(x, y, **kwargs):
    df = kwargs.pop('data')
    y2 = kwargs.pop('y2')
    ax = plt.gca()
    df.plot(x=x, y=y, ax=ax)
    # df.plot(x=x, y=y2, ax=ax)
    # ax.set_ylabel('rate')
    ax2 = ax.twinx()
    color = 'g'
    df.plot(x=x, y=y2, ax=ax2 , color=color)
    ax2.set_ylabel(y2, color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)

def roc_plot_grid(roc_df, fig_path, fig_name):
    # Plot ROC curves
    hue = "serum dilution"
    antigens = natsorted(roc_df['antigen'].unique())
    sns.set_context("notebook")
    assert not roc_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
    palette = sns.color_palette(n_colors=len(roc_df[hue].unique()))
    print('plotting ROC curves...')
    g = sns.FacetGrid(roc_df, hue=hue, col="antigen", col_order=antigens, col_wrap=3, aspect=1.1,
                      hue_kws={'linestyle': ['-', '--', '-.', ':']})
    # g = (g.map(plt.plot, 'False positive rate', 'True positive rate').add_legend())
    g = (g.map_dataframe(roc_plot, 'False positive rate', 'True positive rate', y2='threshold'))
    for antigen, ax in zip(antigens, g.axes.flat):
        sub_df = roc_df[roc_df['antigen'] == antigen]
        auc = sub_df['AUC'].unique()[0]
        ax.set_title(antigen)
        ax.text(0.6, 0.15, 'AUC={:.3f}'.format(auc), fontsize=12)  # add text
        # ax2 = ax.twinx()
        # sub_df.plot(x='False positive rate', y='threshold', ax=ax2)
    plt.savefig(os.path.join(fig_path, fig_name + '.jpg'),
                             dpi=300, bbox_inches='tight')
    plt.close()

def normalize_od_helper(norm_antigen):
    def normalize(df):
        norm_antigen_df = slice_df(df, 'keep', 'antigen', [norm_antigen])
        norm_factor = norm_antigen_df['OD'].mean()
        df['OD'] = df['OD'] / norm_factor
        return df
    return normalize

def normalize_od(df, norm_antigen):
    """fit model to x, y data in dataframe.
    Return a dataframe with fit x, y for plotting
    """
    norm_antigen_df = slice_df(df, 'keep', 'antigen', [norm_antigen])
    df.loc[df['antigen'] == norm_antigen, 'OD'] = norm_antigen_df['OD'] / norm_antigen_df['OD'].mean()
    norm_fn = normalize_od_helper(norm_antigen)
    df = df.groupby(['plate_id']).apply(norm_fn)
    df = df.reset_index()
    return df

def scatter_plot(df,
                 x_col,
                 y_col,
                 title,
                 output_path,
                 output_fname,
                 xlim=None,
                 ylim=None,
                 alpha=1):
    diff_df = df[y_col] - df[x_col]
    me = diff_df.mean()
    mae = diff_df.abs().mean()
    fig = plt.figure()
    fig.set_size_inches((6, 6))
    ax = sns.scatterplot(x=x_col, y=y_col, data=df, alpha=alpha)
    plt.title(title)
    if xlim is None:
        xlim = ax.get_xlim()
    if ylim is None:
        ylim = ax.get_xlim()
    ax.set_xlim(left=xlim[0], right=xlim[1])
    ax.set_ylim(bottom=ylim[0], top=ylim[1])
    xfit = np.linspace(xlim[0], xlim[1], 2)
    plt.plot(xfit, xfit, linewidth=5, color='k', linestyle='--', alpha=0.5)
    ax.text(0.7 * xlim[1], 0.15 * ylim[1], 'Bias={:.3f}'.format(me), fontsize=16)  # add text
    ax.text(0.7 * xlim[1], 0.1 * ylim[1], 'Noise={:.3f}'.format(mae), fontsize=16)  # add text
    plt.savefig(os.path.join(output_path, ''.join([output_fname, '.jpg'])),
                dpi=300, bbox_inches='tight')
    plt.close()

#%% plate 3
# scn_psr_dirs = [r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200630_1647',
#                 r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200630_1657',
#                 r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200630_1708',
#                 r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200630_1739']
# fig_path = os.path.join(r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs',
#                         'pysero_plots')
# scn_psr_plate_ids = ['plate_3'] * len(scn_psr_dirs)
# scn_psr_slice_actions = ['keep', 'drop', 'keep', 'drop']
# scn_psr_well_ids = [['E3', 'E4'], ['E3', 'E4'], ['F2'], ['F2']]
# sera_fit_list = ['Pool']
# sera_cat_list = ['Pool', 'Blank']
# sera_roc_list = sera_cat_list

#%% plate 9
# scn_psr_dirs = [r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-25-28-COVID_June24_OJassay_plate9_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200701_0933',
#                 r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-25-28-COVID_June24_OJassay_plate9_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_1002',
#                 r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-25-28-COVID_June24_OJassay_plate9_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_1028']
# fig_path = os.path.join(r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-25-28-COVID_June24_OJassay_plate9_images/Stitched data from multiple pysero outputs',
#                         'pysero_plots')
# scn_psr_plate_ids = ['plate_9'] * len(scn_psr_dirs)
# scn_psr_slice_actions = [None, 'drop', 'keep']
# scn_psr_well_ids = [None, ['B7', 'B8', 'B10','B11','D3','D8','D9','D10','D11','F8','F9','F11','H2','H7','H8','H9','H10','H11'],
#                 ['B7', 'B8', 'B10','B11','D3','D8','D9','D10','D11','F8','F9','F11','H2','H7','H8','H9','H10','H11']]
# sera_fit_list = ['mab']
# sera_cat_list = ['mab', 'Blank']
# sera_roc_list = sera_cat_list
#%% plate 3 & 9
scn_psr_dirs = [r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200630_1647',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200630_1657',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200630_1708',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200630_1739',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-25-28-COVID_June24_OJassay_plate9_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200701_0933',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-25-28-COVID_June24_OJassay_plate9_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_1002',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-25-28-COVID_June24_OJassay_plate9_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_1028',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-40-02-COVID_June5_OJassay_plate7_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200611_1257',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-40-02-COVID_June5_OJassay_plate7_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200613_1341',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-40-02-COVID_June5_OJassay_plate7_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200709_1024',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-44-32-COVID_June5_OJassay_plate8_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200613_2017',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-44-32-COVID_June5_OJassay_plate8_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200613_2042',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-44-32-COVID_June5_OJassay_plate8_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200613_2051',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-44-32-COVID_June5_OJassay_plate8_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200618_1449',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-21-32-COVID_June24_OJassay_plate4_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200630_1756',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-21-32-COVID_June24_OJassay_plate4_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_0849',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-21-32-COVID_June24_OJassay_plate4_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_0903',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-29-42-COVID_June24_OJassay_plate10_images/Stitched data from multiple pysero outputs/pysero_biotin_fiducial_20200701_1047',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-29-42-COVID_June24_OJassay_plate10_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_1110',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-29-42-COVID_June24_OJassay_plate10_images/Stitched data from multiple pysero outputs/pysero_igg_fiducial_20200701_1141',
                ]

scn_scn_dirs = [r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-18-08-COVID_June24_OJassay_plate3_images/',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-25-28-COVID_June24_OJassay_plate9_images/',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-40-02-COVID_June5_OJassay_plate7_images/',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-05-15-44-32-COVID_June5_OJassay_plate8_images/',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-21-32-COVID_June24_OJassay_plate4_images/',
                r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/2020-06-24-17-29-42-COVID_June24_OJassay_plate10_images/',
                ]

fig_path = os.path.join(r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_scienion/OJ_plate3_9_7_8_4_10',
                        'pysero_plots')
scn_psr_plate_ids = ['plate_3'] * 4 + ['plate_9'] * 3 + ['plate_7'] * 3 + ['plate_8'] * 4 + ['plate_4'] * 3 + ['plate_10'] * 3
scn_scn_plate_ids = ['plate_3', 'plate_9', 'plate_7', 'plate_8', 'plate_4', 'plate_10']
scn_psr_slice_actions = ['keep', 'drop', 'keep', 'drop'] + \
                        [None, 'drop', 'keep'] + \
                        [None, 'drop', 'keep'] + \
                        [None, 'keep', 'drop', 'keep'] + \
                        [None, 'drop', 'keep'] + \
                        [None, 'drop', 'keep']
scn_psr_well_ids = [['E3', 'E4'], ['E3', 'E4'], ['F2'], ['F2']] + \
                   [None, ['B7', 'B8', 'B10','B11','D3','D8','D9','D10','D11','F8','F9','F11','H2','H7','H8','H9','H10','H11'],
                ['B7', 'B8', 'B10','B11','D3','D8','D9','D10','D11','F8','F9','F11','H2','H7','H8','H9','H10','H11']] + \
                   [None, ['B4','F2'], ['B4','F2']] + \
                   [None, ['B1','D4','D6','F5'], ['B1', 'B12','D4','D6','F5','H11'], ['B12','H11']] + \
                   [None, ['B3','B5','D1','D2','D5','F5'], ['B3','B5','D1','D2','D5','F5']] + \
                   [None, ['B12','D2','D3','D1','D5','H2','H5'], ['B12','D1','D5','D2','D3','H2']]
sera_fit_list = ['Pool', 'mab', 'CR3022']
sera_cat_list = ['Pool', 'mab', 'Blank', 'CR3022']
sera_roc_list = sera_cat_list
#%% Nautilus

ntl_dirs = [r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate10_images_655_2020-06-24 18-37-39.863642/0 renamed 16bit/biotin_fiducial/pysero_biotin_fiducial_20200708_1813',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate10_images_655_2020-06-24 18-37-39.863642/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200708_1853',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate10_images_655_2020-06-24 18-37-39.863642/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200709_1431',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate10_images_655_2020-06-24 18-37-39.863642/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200709_1613',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate10_images_655_2020-06-24 18-37-39.863642/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200709_1516',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate3_images_655_2020-06-24 17-58-27.941906/0 renamed 16bit/biotin_fiducial/pysero_biotin_fiducial_20200702_1350',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate3_images_655_2020-06-24 17-58-27.941906/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200727_1954',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate3_images_655_2020-06-24 17-58-27.941906/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200706_0654',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate3_images_655_2020-06-24 17-58-27.941906/0 renamed 16bit/igg_fiducial/pysero_igg_fiducial_20200728_1129',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate9_images_655_2020-06-24 18-26-58.173049/0 renamed/Stitched data/pysero_biotin_fiducial_20200706_0847',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate9_images_655_2020-06-24 18-26-58.173049/0 renamed/Stitched data/pysero_igg_fiducial_20200706_0903',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate9_images_655_2020-06-24 18-26-58.173049/0 renamed/Stitched data/pysero_igg_fiducial_20200706_0937',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate9_images_655_2020-06-24 18-26-58.173049/0 renamed/Stitched data/pysero_igg_fiducial_20200708_1136',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate7_images_2020-06-08 19-54-45.277504/0 renamed/biotin fiducial/pysero_biotin fiducial_20200727_2016',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate7_images_2020-06-08 19-54-45.277504/0 renamed/igg fiducial/pysero_igg fiducial_20200727_2038',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate7_images_2020-06-08 19-54-45.277504/0 renamed/igg fiducial/pysero_igg fiducial_20200728_0822',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate7_images_2020-06-08 19-54-45.277504/0 renamed/igg fiducial/pysero_igg fiducial_20200728_0902',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate8_images_2020-06-08 20-38-40.554230/0 renamed/biotin_fiducial/pysero_biotin_fiducial_20200728_1018',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate8_images_2020-06-08 20-38-40.554230/0 renamed/igg_fiducial/pysero_igg_fiducial_20200728_1044',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-08-COVID_June5_OJassay_plate8_images_2020-06-08 20-38-40.554230/0 renamed/igg_fiducial/pysero_igg_fiducial_20200728_1108',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate4_images_655_2020-06-24 18-13-40.894052/0 renamed/biotin_fiducial/pysero_biotin_fiducial_20200728_0921',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate4_images_655_2020-06-24 18-13-40.894052/0 renamed/igg_fiducial/pysero_igg_fiducial_20200709_1411',
            r'/Volumes/GoogleDrive/My Drive/ELISAarrayReader/images_nautilus/2020-06-24-COVID_June24_OJAssay_Plate4_images_655_2020-06-24 18-13-40.894052/0 renamed/igg_fiducial/pysero_igg_fiducial_20200728_0956',]

ntl_plate_ids = ['plate_10'] * 5 + ['plate_3'] * 4 + ['plate_9'] * 4 + ['plate_7'] * 4 + ['plate_8'] * 3 + ['plate_4'] * 3

ntl_slice_actions = [None, 'drop','drop','keep','keep'] +\
                  [None,'drop','keep','keep'] +\
                  [None, 'drop','drop','keep'] +\
                  [None, 'keep', 'drop', 'keep'] +\
                  [None,'drop','keep'] +\
                   [None,'drop','keep']


ntl_well_ids = [None,['B2', 'B6', 'B7', 'B8', 'B9', 'B10', 'B11', 'B12', 'D6', 'D7', 'D8', 'D9', 'D10', 'D11', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11'],['B1', 'B2', 'B4', 'B12', 'D2', 'D4', 'D5', 'D8', 'F3','H4','B3','B5','D1','D2','D3','D12','F1','F2','F4','F5','F12','H1','H2','H3','H12'],['B12', 'H4'],['B2','D8']] + \
           [None, ['D4','D7','D8','D9','D10', 'D11', 'F7', 'F8', 'F9','F10', 'F12','H7','H9','H10','H11'],['D8','D9','D10','D11','F8','F9','F10','H7','H9','H10','H11'], ['D4','D7','F7','F12']] + \
            [None, ['B1',' B2', 'B5','B7','B8','b9','b10','B11','D3','D5','D7','D8','D9','D10','D11','F7','F8','F9','F10','F11','H7','H8','H9','H10','H11'], ['B5','D5','D6'],['B5','D2','D6']] + \
            [None, ['B2','B12','D3','D12', 'F4','F12','H1','H2','H3','H6'], ['B2','B12','D3','D12', 'F4','F12','H1','H2','H3','H6','B7','B8','B9','B10','B11','D7','D8','D9','D10','D11','F2','F6','F7','F8','F9','F10','F11','H4','H7','H8','H9','H10','H11','H12'], ['B7','B8','B9','B10','B11','D7','D8','D9','D10','D11','F2','F6','F7','F8','F9', 'F10','F11','H4','H7','H8','H9','H10','H11','H12']] + \
            [None, ['B9','B10','B11','B12','D7','D8','D9','D11','F6','F7','F8','F9','F10','F11','H6','H7','H8','H9', 'H11','H12'], ['B9','B10','B11','D7','D8','D9','D11','F6','F7','F8','F9','F10','F11','H6','H7','H8','H9', 'H11','H12']] + \
            [None,['D1','D5','D12','H2'], ['D1','D5','D12','H2']]
#%%
os.makedirs(fig_path, exist_ok=True)
df_list = []
scn_df = pd.DataFrame()

for scn_scn_dir, plate_id, in zip(scn_scn_dirs, scn_scn_plate_ids):
    metadata_path = os.path.join(scn_scn_dir, 'pysero_output_data_metadata.xlsx')
    with pd.ExcelFile(metadata_path) as meta_file:
        antigen_df = read_antigen_info(meta_file)
        plate_info_df = read_plate_info(meta_file)
    plate_info_df['plate_id'] = plate_id
    scn_fname = [f for f in os.listdir(scn_scn_dir) if '_analysis.xlsx' in f]
    scn_path = os.path.join(scn_scn_dir, scn_fname[0])
    scn_df_tmp = read_scn_output(scn_path, plate_info_df)
    # Join Scienion data with plateInfo
    scn_df_tmp = pd.merge(scn_df_tmp,
                      antigen_df,
                      how='left', on=['antigen_row', 'antigen_col'])
    scn_df_tmp = pd.merge(scn_df_tmp,
                      plate_info_df,
                      how='right', on=['well_id'])
    scn_df = scn_df.append(scn_df_tmp, ignore_index=True)
scn_df['pipeline'] = 'scienion'
scn_df.dropna(subset=['OD'], inplace=True)
#%%
for data_folder, slice_action, well_id, plate_id in \
        zip(scn_psr_dirs, scn_psr_slice_actions, scn_psr_well_ids, scn_psr_plate_ids):
    metadata_path = os.path.join(data_folder, 'pysero_output_data_metadata.xlsx')
    OD_path = os.path.join(data_folder, 'median_ODs.xlsx')
    int_path = os.path.join(data_folder, 'median_intensities.xlsx')
    bg_path = os.path.join(data_folder, 'median_backgrounds.xlsx')

    with pd.ExcelFile(metadata_path) as meta_file:
        antigen_df = read_antigen_info(meta_file)
        plate_info_df = read_plate_info(meta_file)
    plate_info_df['plate_id'] = plate_id
    OD_df = read_pysero_output(OD_path, antigen_df, file_type='od')
    int_df = read_pysero_output(int_path, antigen_df, file_type='int')
    bg_df = read_pysero_output(bg_path, antigen_df, file_type='bg')
    OD_df = pd.merge(OD_df,
                     antigen_df[['antigen_row', 'antigen_col', 'antigen type']],
                     how='left', on=['antigen_row', 'antigen_col'])
    OD_df = pd.merge(OD_df,
                     plate_info_df,
                     how='right', on=['well_id'])
    pysero_df = pd.merge(OD_df,
                         int_df,
                         how='left', on=['antigen_row', 'antigen_col', 'well_id'])
    pysero_df = pd.merge(pysero_df,
                         bg_df,
                         how='left', on=['antigen_row', 'antigen_col', 'well_id'])
    pysero_df['pipeline'] = 'python'
    # pysero_df = pysero_df.append(scn_df)
    pysero_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    pysero_df.dropna(subset=['OD'], inplace=True)
    pysero_df = slice_df(pysero_df, slice_action, 'well_id', well_id)
    df_list.append(pysero_df)
#%%
for data_folder, slice_action, well_id, plate_id in \
        zip(ntl_dirs, ntl_slice_actions, ntl_well_ids, ntl_plate_ids):
    metadata_path = os.path.join(data_folder, 'pysero_output_data_metadata.xlsx')
    OD_path = os.path.join(data_folder, 'median_ODs.xlsx')
    int_path = os.path.join(data_folder, 'median_intensities.xlsx')
    bg_path = os.path.join(data_folder, 'median_backgrounds.xlsx')

    with pd.ExcelFile(metadata_path) as meta_file:
        antigen_df = read_antigen_info(meta_file)
        plate_info_df = read_plate_info(meta_file)
    plate_info_df['plate_id'] = plate_id
    OD_df = read_pysero_output(OD_path, antigen_df, file_type='od')
    int_df = read_pysero_output(int_path, antigen_df, file_type='int')
    bg_df = read_pysero_output(bg_path, antigen_df, file_type='bg')
    OD_df = pd.merge(OD_df,
                     antigen_df[['antigen_row', 'antigen_col', 'antigen type']],
                     how='left', on=['antigen_row', 'antigen_col'])
    OD_df = pd.merge(OD_df,
                     plate_info_df,
                     how='right', on=['well_id'])
    pysero_df = pd.merge(OD_df,
                         int_df,
                         how='left', on=['antigen_row', 'antigen_col', 'well_id'])
    pysero_df = pd.merge(pysero_df,
                         bg_df,
                         how='left', on=['antigen_row', 'antigen_col', 'well_id'])
    pysero_df['pipeline'] = 'nautilus'
    pysero_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    pysero_df.dropna(subset=['OD'], inplace=True)
    pysero_df = slice_df(pysero_df, slice_action, 'well_id', well_id)
    df_list.append(pysero_df)
df_list.append(scn_df)
#%% Concatenate dataframes
stitched_pysero_df = pd.concat(df_list)
stitched_pysero_df.reset_index(drop=True, inplace=True)
# remove empty xkappa-biotin spots, round off dilution
stitched_pysero_df = stitched_pysero_df[(stitched_pysero_df['antigen'] != 'xkappa-biotin') |
                        (stitched_pysero_df['antigen type'] == 'Fiducial')]
stitched_pysero_df['serum dilution'] = stitched_pysero_df['serum dilution'].round(7)
#%% 4PL fit
slice_cols = ['pipeline', 'serum ID']
slice_keys = [['python'], sera_fit_list]
scn_psr_slice_actions = ['keep', 'keep']
serum_df = stitched_pysero_df.copy()
serum_df_fit = stitched_pysero_df.copy()
for col, action, key in zip(slice_cols, scn_psr_slice_actions, slice_keys):
    serum_df = slice_df(serum_df, action, col, key)
    serum_df_fit = slice_df(serum_df_fit, action, col, key)
serum_df_fit = fit2df(serum_df_fit, fourPL, 'python')

#%%
def roc_plot(x, y, **kwargs):
    df = kwargs.pop('data')
    y2 = kwargs.pop('y2')
    ax = plt.gca()
    df.plot(x=x, y=y, ax=ax, legend=False)
    # df.plot(x=x, y=y2, ax=ax)
    # ax.set_ylabel('rate')
    ax2 = ax.twinx()
    color = 'g'
    df.plot(x=x, y=y2, ax=ax2, color=color, style='--', legend=False)
    ax2.set_ylabel(y2, color=color)  # we already handled the x-label with ax1
    ax2.tick_params(axis='y', labelcolor=color)

def roc_plot_grid(roc_df, fig_path, fig_name):
    # Plot ROC curves
    hue = "serum dilution"
    antigens = natsorted(roc_df['antigen'].unique())
    sns.set_context("notebook")
    assert not roc_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
    palette = sns.color_palette(n_colors=len(roc_df[hue].unique()))
    print('plotting ROC curves...')
    g = sns.FacetGrid(roc_df, hue=hue, col="antigen", col_order=antigens, col_wrap=3, aspect=1.2,
                      hue_kws={'linestyle': ['-', '--', '-.', ':']})
    # g = (g.map(plt.plot, 'False positive rate', 'True positive rate').add_legend())
    g = (g.map_dataframe(roc_plot, 'False positive rate', 'True positive rate', y2='threshold'))
    for antigen, ax in zip(antigens, g.axes.flat):
        sub_df = roc_df[roc_df['antigen'] == antigen]
        auc = sub_df['AUC'].unique()[0]
        ax.set_title(antigen)
        ax.text(0.6, 0.15, 'AUC={:.3f}'.format(auc), fontsize=12)  # add text
        # ax2 = ax.twinx()
        # sub_df.plot(x='False positive rate', y='threshold', ax=ax2)
    plt.savefig(os.path.join(fig_path, fig_name + '.jpg'),
                             dpi=300, bbox_inches='tight')
    plt.close()

def thr_plot_grid(roc_df, fig_path, fig_name):
    # Plot ROC curves
    hue = 'category'
    antigens = natsorted(roc_df['antigen'].unique())
    sns.set_context("notebook")
    assert not roc_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
    palette = sns.color_palette(n_colors=len(roc_df[hue].unique()))
    print('plotting ROC curves...')
    g = sns.FacetGrid(roc_df, hue=hue, col="antigen", col_order=antigens, col_wrap=3, aspect=1.05,
                      hue_kws={'linestyle': ['-', '--', '-.', ':']})
    g = (g.map(plt.plot, 'threshold', 'rate').add_legend())
    # g = (g.map_dataframe(roc_plot, 'False positive rate', 'True positive rate', y2='threshold'))
    for antigen, ax in zip(antigens, g.axes.flat):
        sub_df = roc_df[roc_df['antigen'] == antigen]
        auc = sub_df['AUC'].unique()[0]
        ax.set_title(antigen)
        ax.text(1.5, 0.2, 'AUC={:.3f}'.format(auc), fontsize=12)  # add text
        # ax2 = ax.twinx()
        # sub_df.plot(x='False positive rate', y='threshold', ax=ax2)
    plt.savefig(os.path.join(fig_path, fig_name + '.jpg'),
                             dpi=300, bbox_inches='tight')
    plt.close()
#%% functions to compute ROC curves and AUC
# for plate_id in stitched_pysero_df['plate_id'].unique():
# for plate_id in ['plate_8']:
#     slice_cols = ['pipeline', 'serum ID', 'plate_id']
#     slice_keys = [['python'], sera_roc_list, [plate_id]]
#     scn_psr_slice_actions = ['keep', 'drop', 'keep']
for pipeline in stitched_pysero_df['pipeline'].unique():
    slice_cols = ['pipeline', 'serum ID']
    slice_keys = [[pipeline], sera_roc_list]
    scn_psr_slice_actions = ['keep', 'drop']
    roc_df = stitched_pysero_df.copy()
    for col, action, key in zip(slice_cols, scn_psr_slice_actions, slice_keys):
        roc_df = slice_df(roc_df, action, col, key)
    roc_df = get_roc_df(roc_df)
    # roc_plot_grid(roc_df, fig_path, 'ROC')
    # roc_plot_grid(roc_df, fig_path, 'ROC_' + plate_id)
    roc_plot_grid(roc_df, fig_path, 'ROC_' + pipeline)
    roc_df = roc_df.melt(id_vars=['antigen',
                                 'secondary ID',
                                 'secondary dilution',
                                 'serum dilution',
                                 'pipeline',
                                'threshold',
                                'AUC'],
                         var_name='category',
                         value_name='rate'
                         )
    # thr_plot_grid(roc_df, fig_path, 'ROC_thr')
    # thr_plot_grid(roc_df, fig_path, 'ROC_thr_' + plate_id)
    thr_plot_grid(roc_df, fig_path, 'ROC_thr_' + pipeline)
#%% Plot categorical scatter plot for episurvey
for plate_id in stitched_pysero_df['plate_id'].unique():
# for plate_id in ['plate_4']:
    slice_cols = ['pipeline', 'serum ID', 'plate_id']
    slice_keys = [['python'], sera_cat_list, [plate_id]]
    scn_psr_slice_actions = ['keep', 'drop', 'keep']
    antigens = natsorted(stitched_pysero_df['antigen'].unique())
    serum_df = stitched_pysero_df.copy()
    for col, action, key in zip(slice_cols, scn_psr_slice_actions, slice_keys):
        serum_df = slice_df(serum_df, action, col, key)
    assert not serum_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
    # Draw a categorical scatterplot to show each observation
    g = sns.catplot(x="serum type", y="OD", hue="serum type", col_order=antigens, col="antigen",kind="swarm",
                    palette=["r", "c", "y"], data=serum_df, col_wrap=5)
    plt.savefig(os.path.join(fig_path, 'catplot_{}.jpg'.format(plate_id)),
                                  dpi=300, bbox_inches='tight')
    g.set(ylim=(-0.05, 0.4))
    plt.savefig(os.path.join(fig_path, 'catplot_zoom_{}.jpg'.format(plate_id)),
                                  dpi=300, bbox_inches='tight')
    plt.close('all')
#%% pivot the dataframe for xy scatter plot
norm_antigen = 'xIgG Fc'
# norm_antigen = 'xkappa-biotin'
suffix = '_'.join([norm_antigen, 'norm_per_plate'])
# norm_antigen = ''
df_norm = stitched_pysero_df.copy()
df_norm = normalize_od(stitched_pysero_df, norm_antigen)
pysero_df_pivot = pd.pivot_table(df_norm, values='OD',
                             index=['well_id', 'antigen_row', 'antigen_col', 'serum ID', 'secondary ID', 'secondary dilution',
       'serum type', 'serum dilution', 'antigen', 'antigen type', 'pipeline'],
                             columns=['plate_id'])
pysero_df_pivot.reset_index(inplace=True)
#%%
limit = [0, 2.0]
x_cols = ['plate_3', 'plate_7', 'plate_4']
y_cols = ['plate_9', 'plate_8', 'plate_10']

# suffix = ''
antigen_OD_df = slice_df(pysero_df_pivot, 'keep', 'antigen type', ['Diagnostic'])
biotin_OD_df = slice_df(pysero_df_pivot, 'keep', 'antigen', ['xkappa-biotin'])
# biotin_OD_df = slice_df(biotin_OD_df, 'keep', 'antigen type', ['Fiducial'])
igg_OD_df = slice_df(pysero_df_pivot, 'keep', 'antigen', ['xIgG Fc'])
#%%
for x_col, y_col in zip(x_cols, y_cols):
    # scatter_plot(pysero_df_pivot, x_col, y_col, fig_path, '_'.join(['OD_scatter', x_col, y_col, suffix]), xlim=limit, ylim=limit)
    scatter_plot(antigen_OD_df, x_col, y_col, 'antigen', fig_path,
                 '_'.join(['antigen_OD_scatter', x_col, y_col, suffix]), xlim=limit, ylim=limit)
    scatter_plot(biotin_OD_df, x_col, y_col, 'biotin', fig_path,
                 '_'.join(['biotin_OD_scatter', x_col, y_col, suffix]), xlim=limit, ylim=limit)
    scatter_plot(igg_OD_df, x_col, y_col, 'igg', fig_path,
                 '_'.join(['igg_OD_scatter', x_col, y_col, suffix]), xlim=limit, ylim=limit)
#%% functions to compute ROC curves and AUC
slice_cols = ['pipeline', 'serum ID']
slice_keys = [['python'], sera_roc_list]
scn_psr_slice_actions = ['keep', 'drop']
roc_df = df_norm.copy()
for col, action, key in zip(slice_cols, scn_psr_slice_actions, slice_keys):
    roc_df = slice_df(roc_df, action, col, key)
roc_df = get_roc_df(roc_df)
roc_plot_grid(roc_df, fig_path, '_'.join(['ROC', suffix]))
#%%
serum_df = serum_df[(serum_df['antigen'] == 'SARS CoV2 RBD 500')]
neg_max = serum_df[(serum_df['serum type'] == 'negative')]['OD'].max()
overlap_df = serum_df[(serum_df['serum type'] == 'positive') & (serum_df['OD'] <= neg_max)]
#%%
serum_id = stitched_pysero_df['serum ID'].unique()
#%% plot the ODs and fits
sera_4pl_list = [' '.join([x, 'fit']) for x in sera_fit_list]
markers = 'o'
hue = 'serum ID'
style = 'serum type'
antigens = natsorted(stitched_pysero_df['antigen'].unique())
assert not serum_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
palette = sns.color_palette(n_colors=len(serum_df[hue].unique()))
print('plotting standard curves...')
g = sns.lmplot(x="serum dilution", y="OD", col_order=antigens,
                hue=hue, hue_order=sera_fit_list, col="antigen", ci='sd', palette=palette, markers=markers,
                 data=serum_df, col_wrap=5, fit_reg=False, x_estimator=np.mean)
palette = sns.color_palette(n_colors=len(serum_df_fit[hue].unique()))
for antigen, ax in zip(antigens, g.axes.flat):
    df_fit = serum_df_fit[(serum_df_fit['antigen'] == antigen)]
    sns.lineplot(x="serum dilution", y="OD", hue=hue, hue_order=sera_4pl_list, data=df_fit,
                 style=style, palette=palette,
                 ax=ax, legend=False)
    ax.set(xscale="log")
plt.savefig(os.path.join(fig_path, '{}_fit.jpg'.format(sera_fit_list)),dpi=300, bbox_inches='tight')

for antigen, ax in zip(antigens, g.axes.flat):
    ax.set(ylim=[-0.05, 1.5])
plt.savefig(os.path.join(fig_path, '{}_fit_zoom.jpg'.format(sera_fit_list)),dpi=300, bbox_inches='tight')
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
# python_sera_df = pysero_df[pysero_df['serum ID'].isin(sera_list) & (pysero_df['pipeline']==pipeline)]
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
# for serum in pysero_df['serum ID'].unique():
#     serum_df = pysero_df[pysero_df['serum ID']==serum]
#     for y in y_list:
#         g = sns.relplot(x="serum dilution", y=y,
#                     hue="pipeline", col="antigen", ci='sd', style='pipeline',
#                      data=serum_df, col_wrap=3, estimator=np.mean, kind='line')
#         g.set(xscale="log", ylim=[-0.05, 1.2])
#         plt.savefig(os.path.join(fig_path, '_'.join([serum, y + '.jpg'])), dpi=300, bbox_inches='tight')
#         plt.close()
#
#





