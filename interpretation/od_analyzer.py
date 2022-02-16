import pandas as pd
import numpy as np
import os
import re
from matplotlib import pyplot as plt
import seaborn as sns
from interpretation.plotting import roc_plot_grid, total_plots
from interpretation.report_reader import slice_df, normalize_od, read_output_batch
import array_analyzer.extract.constants as constants


def read_config(input_dir):
    """
    Load analysis config from the input directory.
    :param str input_dir: input directory
    :return dataframe ntl_dirs_df: 'multisero output dirs' tab in the config file.
    :return dataframe scn_scn_df: 'scienion output dirs' tab in the config file.
    :return dataframe plot_setting_df: 'general plotting settings' tab in the config file.
    :return dataframe roc_param_df: 'ROC plot' tab in the config file.
    :return dataframe cat_param_df: 'categorical plot' tab in the config file.
    :return dataframe fit_param_df: 'standard curves' tab in the config file.
    """
    config = dict()
    config_path = os.path.join(input_dir, constants.METADATA_FILE)
    if constants.METADATA_FILE not in os.listdir(input_dir):
        raise IOError("analysis config file not found, aborting")
    # check that the xlsx file contains necessary worksheets
    ntl_dirs_df = scn_scn_df = roc_param_df = cat_param_df = fit_param_df = pd.DataFrame()
    with pd.ExcelFile(config_path) as config_file:
        plot_setting_df = pd.read_excel(config_file, sheet_name='general plotting settings',
                                        index_col=0, squeeze=True, usecols='A,B')
        plot_setting_df.where(plot_setting_df.notnull(), None, inplace=True)
        if 'ROC plot' in config_file.sheet_names:
            roc_param_df = pd.read_excel(config_file, sheet_name='ROC plot',
                                         index_col=0, squeeze=True, usecols='A,B')
            # replace NaN with None
            roc_param_df.where(roc_param_df.notnull(), None, inplace=True)
            if roc_param_df['serum ID'] is not None:
                roc_param_df['serum ID'] = re.split(r'\s*,\s*', roc_param_df['serum ID'])
        if 'categorical plot' in config_file.sheet_names:
            cat_param_df = pd.read_excel(config_file, sheet_name='categorical plot',
                                         index_col=0, squeeze=True, usecols='A,B')
            cat_param_df.where(cat_param_df.notnull(), None, inplace=True)
            if cat_param_df['serum ID'] is not None:
                cat_param_df['serum ID'] = re.split(r'\s*,\s*', cat_param_df['serum ID'])
        if 'standard curves' in config_file.sheet_names:
            fit_param_df = pd.read_excel(config_file, sheet_name='standard curves',
                                         index_col=0, squeeze=True, usecols='A,B')
            fit_param_df.where(fit_param_df.notnull(), None, inplace=True)
            if fit_param_df['serum ID'] is not None:
                fit_param_df['serum ID'] = re.split(r'\s*,\s*', fit_param_df['serum ID'])
        if not constants.LOAD_REPORT:
            assert ('multisero output dirs' in config_file.sheet_names) or \
            ('scienion output dirs' in config_file.sheet_names), \
            "sheet by name 'multisero output dirs' or 'scienion output dirs' are required " \
            "in analysis config file when load_report is False, aborting"
            if 'multisero output dirs' in config_file.sheet_names:
                ntl_dirs_df = pd.read_excel(config_file, sheet_name='multisero output dirs', comment='#')
                if not ntl_dirs_df.isna().loc[0, 'well ID']:
                    ntl_dirs_df['well ID'] = ntl_dirs_df['well ID'].str.split(pat=r'\s*,\s*')
            if 'scienion output dirs' in config_file.sheet_names:
                scn_scn_df = pd.read_excel(config_file, sheet_name='scienion output dirs', comment='#')
    return ntl_dirs_df, scn_scn_df, plot_setting_df, roc_param_df, cat_param_df, fit_param_df

def analyze_od(input_dir, output_dir, load_report):
    """
    Perform analysis on multisero or scienion OD outputs specified in the config files.
    Save the combined table as 'master report' in the output directory.
    :param str input_dir: Input directory
    :param str output_dir: Output directory
    :param bool load_report: If True, load the saved 'master report' in the output directory
    from the previous run. Load from the master report is much faster.
    """
    os.makedirs(output_dir, exist_ok=True)
    ntl_dirs_df, scn_scn_df, plot_setting_df, roc_param_df, cat_param_df, fit_param_df =\
        read_config(input_dir)
    stitched_multisero_df = read_output_batch(output_dir, ntl_dirs_df, scn_scn_df, load_report)
    if plot_setting_df['antigens to plot'] == 'all':
        plot_setting_df['antigens to plot'] = stitched_multisero_df['antigen'].unique()
    split_plots_by = plot_setting_df['split plots by']
    split_plots_vals = [None]
    if split_plots_by is not None:
        split_plots_vals = stitched_multisero_df[split_plots_by].unique()
    norm_antigen = plot_setting_df['normalize OD by']
    norm_group = 'plate'
    aggregate = 'mean'
    # aggregate = None
    #antigen_list = plot_setting_df['antigens to plot']
    antigen_list = []
    #antigen_list.append(plot_setting_df['antigens to plot'])
    antigen_list.append('xIgG Fc')
    antigen_list.append('xIgG Fc Old D12')
    #antigen_list.append('DENV1 100 RVP')
    antigen_list.append('DENV1 125 NS1 Old D12')
    antigen_list.append('DENV1 125 NS1 New D12')
    #antigen_list.append('DENV1 100 EDIII')
    #antigen_list.append('DENV2 50 VLP')
    #antigen_list.append('DENV2 100 RVP')
    antigen_list.append('DENV2 250 NS1 Old D12')
    antigen_list.append('DENV2 250 NS1 New D12')
    #antigen_list.append('DENV2 100 EDIII')
    #antigen_list.append('DENV3 50 VLP')
    #antigen_list.append('DENV3 100 RVP')
    #antigen_list.append('DENV3 250 NS1')
    #antigen_list.append('DENV3 100 EDIII')
    #antigen_list.append('xmouse IgG Fc')
    #antigen_list.append('ZIKA 225 NS1')
    #antigen_list.append('CHIK max VLP')
    #antigen_list.append('JEV 245 NS1')
    #antigen_list.append('DENV4 100 EDIII')
    #antigen_list.append('DENV2 50 VLP')
    suffix = ''
    df_norm = normalize_od(stitched_multisero_df.copy(), norm_antigen, group=norm_group)
    if norm_antigen is not None:
        suffix = '_'.join([suffix, norm_antigen, 'norm_per', norm_group])

    if aggregate is not None:
        df_norm = df_norm.groupby(['antigen', 'antigen type', 'serum ID', 'well_id', 'plate ID', 'sample type',
                                 'serum type', 'serum dilution', 'serum cat', 'pipeline', 'secondary ID',
                                 'secondary dilution','PRNT'])['OD'].mean().reset_index()
        suffix = '_'.join([suffix, aggregate])

    for split_val in split_plots_vals:
        split_suffix = suffix
        if split_val is not None:
            split_suffix = '_'.join([split_suffix, split_val])
        df_norm_sub = slice_df(df_norm, 'keep', split_plots_by, split_val)
        slice_cols = [split_plots_by, 'antigen type', 'antigen']
        slice_keys = [[split_val], ['Diagnostic','Positive'], antigen_list]
        #slice_keys = [[split_val], ['Negative'], antigen_list]
        slice_actions = ['keep', 'keep', 'keep']
        # general slicing
        for col, action, key in zip(slice_cols, slice_actions, slice_keys):
            df_norm_sub = slice_df(df_norm_sub, action, col, key)
        #%% compute ROC curves and AUC
        if not roc_param_df.empty:
            sera_roc_list = roc_param_df['serum ID']
            slice_action = roc_param_df['serum ID action']
            # plot specific slicing
            roc_df = slice_df(df_norm_sub, slice_action, 'serum ID', sera_roc_list)
            fpr = 1 - roc_param_df['specificity']
            ci = roc_param_df['confidence interval']
            hue = roc_param_df['hue']
            # df_norm = offset_od(df_norm, offset_antigen, offset_group)
            roc_suffix = split_suffix
            if ci is not None:
                roc_suffix = '_'.join([roc_suffix, 'ci'])
            #%%
            print('{} unique positive sera'.format(len(roc_df.loc[roc_df['serum type']=='positive', 'serum ID'].unique())))
            print('{} unique negative sera'.format(len(roc_df.loc[roc_df['serum type'] == 'negative', 'serum ID'].unique())))
            roc_plot_grid(roc_df, constants.RUN_PATH, '_'.join(['ROC', roc_suffix]), 'png', ci=ci, fpr=fpr, hue=hue)
    #%% Plot categorical scatter plot for episurvey
        if not cat_param_df.empty:
            sera_cat_list = cat_param_df['serum ID']
            slice_action = cat_param_df['serum ID action']
            split_subplots_by = cat_param_df['split subplots by']
            hue = cat_param_df['hue']
            # plot specific slicing
            cat_df = slice_df(df_norm_sub, slice_action, 'serum ID', sera_cat_list) #serum ID --> antigen
            assert not cat_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
            sns.set_context("talk")
            g = sns.catplot(x="serum type", y="OD", hue=hue, col=split_subplots_by, kind="swarm",
                            data=cat_df, col_wrap=3)
            g.set_xticklabels(rotation=65, horizontalalignment='right')


            plt.savefig(os.path.join(constants.RUN_PATH, 'catplot_{}.png'.format(split_suffix)),
                                          dpi=300, bbox_inches='tight')
            if cat_param_df['zoom']:
                g.set(ylim=(-0.05, 0.4))
                plt.savefig(os.path.join(constants.RUN_PATH, 'catplot_zoom_{}.png'.format(split_suffix)),
                                              dpi=300, bbox_inches='tight')
        #%% 4PL fit
        if not fit_param_df.empty:
            slice_action = fit_param_df['serum ID action']
            hue = fit_param_df['hue']
            dilution_df = slice_df(df_norm_sub, slice_action, 'serum ID', fit_param_df['serum ID'])
            split_subplots_by = fit_param_df['split subplots by']
            total_plots(dilution_df, constants.RUN_PATH, 'fit_{}'.format(split_suffix), 'png', hue=hue,
                                zoom=fit_param_df['zoom'], split_subplots_by=split_subplots_by, col_wrap=2)
        plt.close('all')

