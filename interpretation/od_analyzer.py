import pandas as pd
import numpy as np
import os
import re
import matplotlib
from matplotlib import pyplot as plt
from natsort import natsorted
import seaborn as sns
from interpretation.plotting import fourPL, fit2df, roc_plot_grid, joint_plot
from interpretation.report_reader import read_plate_info, read_antigen_info, read_pysero_output, read_scn_output, \
    slice_df, normalize_od, offset_od
import array_analyzer.extract.constants as constants

def read_config(input_dir):
    """
    Load analysis config from the input directory.
    :param str input_dir: input directory
    :return dataframe :
    """
    config_path = os.path.join(input_dir, constants.METADATA_FILE)
    if constants.METADATA_FILE not in os.listdir(input_dir):
        raise IOError("analysis config file not found, aborting")
    # check that the xlsx file contains necessary worksheets
    ntl_dirs_df = scn_scn_df = roc_param_df = cat_param_df = fit_param_df = pd.DataFrame()
    with pd.ExcelFile(config_path) as config_file:
        plot_setting_df = pd.read_excel(config_file, sheet_name='general plotting settings',
                                        index_col=0, squeeze=True)
        if 'ROC plot' in config_file.sheet_names:
            roc_param_df = pd.read_excel(config_file, sheet_name='ROC plot',
                                         index_col=0, squeeze=True)
            # replace NaN with None
            roc_param_df.where(roc_param_df.notnull(), None, inplace=True)
            roc_param_df['serum ID'] = re.split(r'\s*,\s*', roc_param_df['serum ID'])
        if 'categorical plot' in config_file.sheet_names:
            cat_param_df = pd.read_excel(config_file, sheet_name='categorical plot',
                                         index_col=0, squeeze=True)
            cat_param_df.where(cat_param_df.notnull(), None, inplace=True)
            cat_param_df['serum ID'] = re.split(r'\s*,\s*', cat_param_df['serum ID'])
        if 'standard curves' in config_file.sheet_names:
            fit_param_df = pd.read_excel(config_file, sheet_name='standard curves',
                                         index_col=0, squeeze=True)
            fit_param_df.where(fit_param_df.notnull(), None, inplace=True)
            fit_param_df['serum ID'] = re.split(r'\s*,\s*', fit_param_df['serum ID'])
        if not constants.LOAD_REPORT:
            assert ('pysero output dirs' in config_file.sheet_names) or \
            ('scienion output dirs' in config_file.sheet_names), \
            "sheet by name 'pysero output dirs' or 'scienion output dirs' are required " \
            "in analysis config file when load_report is False, aborting"
            if 'pysero output dirs' in config_file.sheet_names:
                ntl_dirs_df = pd.read_excel(config_file, sheet_name='pysero output dirs')
                ntl_dirs_df['well ID'] = ntl_dirs_df['well ID'].str.split(pat=r'\s*,\s*')
            if 'scienion output dirs' in config_file.sheet_names:
                scn_scn_df = pd.read_excel(config_file, sheet_name='scienion output dirs')
    return ntl_dirs_df, scn_scn_df, plot_setting_df, roc_param_df, cat_param_df, fit_param_df

def read_od_batch(output_dir, ntl_dirs_df, scn_scn_df, load_report):
    if not load_report:
        df_list = []
        scn_df = pd.DataFrame()
        if not scn_scn_df.empty:
            for scn_scn_dir, plate_id, in zip(scn_scn_df['directory'], scn_scn_df['plate ID']):
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

        if not ntl_dirs_df.empty:
            for data_folder, slice_action, well_id, plate_id in \
                    zip(ntl_dirs_df['directory'], ntl_dirs_df['well action'],
                        ntl_dirs_df['well ID'], ntl_dirs_df['plate ID']):
                print('Load {}...'.format(data_folder))
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
        #% Concatenate dataframes
        stitched_pysero_df = pd.concat(df_list)
        stitched_pysero_df.reset_index(drop=True, inplace=True)
        # remove empty xkappa-biotin spots, round off dilution
        stitched_pysero_df = stitched_pysero_df[(stitched_pysero_df['antigen'] != 'xkappa-biotin') |
                                (stitched_pysero_df['antigen type'] == 'Fiducial')]
        stitched_pysero_df['serum dilution'] = stitched_pysero_df['serum dilution'].round(7)
        stitched_pysero_df.to_csv(os.path.join(output_dir, 'master_report.csv'))
    else:
        stitched_pysero_df = pd.read_csv(os.path.join(output_dir, 'master_report.csv'), index_col=0, low_memory=False)
    return stitched_pysero_df

def analyze_od(input_dir, output_dir, load_report):
    os.makedirs(output_dir, exist_ok=True)
    scn_scn_dirs = scn_scn_plate_ids = ntl_dirs = \
        ntl_slice_actions = ntl_well_ids = ntl_plate_ids = None
    ntl_dirs_df, scn_scn_df, plot_setting_df, roc_param_df, cat_param_df, fit_param_df =\
        read_config(input_dir)
    stitched_pysero_df = read_od_batch(output_dir, ntl_dirs_df, scn_scn_df, load_report)
    # fix metadata error
    stitched_pysero_df.loc[stitched_pysero_df['antigen'] == 'xIgG Fc', 'antigen type'] = 'Positive'
    if plot_setting_df['antigens to plot'] == 'all':
        plot_setting_df['antigens to plot'] = stitched_pysero_df['antigen'].unique()
#%%
    split_cols = plot_setting_df['split plots by']
    norm_antigen = plot_setting_df['normalize OD by']
    norm_group = 'plate'
    aggregate = 'mean'
    # aggregate = None
    antigen_list = plot_setting_df['antigens to plot']
    suffix = ''
    df_norm = normalize_od(stitched_pysero_df.copy(), norm_antigen, group=norm_group)
    if norm_antigen is not None:
        suffix = '_'.join([suffix, norm_antigen, 'norm_per', norm_group])

    if aggregate is not None:
        df_norm = df_norm.groupby(['antigen', 'antigen type', 'serum ID', 'well_id', 'plate_id', 'sample type',
                                 'serum type', 'serum dilution', 'pipeline', 'secondary ID',
                                 'secondary dilution'])['OD'].mean().reset_index()
        suffix = '_'.join([suffix, aggregate])

    for split_val in stitched_pysero_df[split_cols].unique():
        suffix = '_'.join([suffix, split_val])
        df_norm_sub = slice_df(df_norm, 'keep', split_cols, [split_val])
        slice_cols = [split_cols, 'antigen type', 'antigen']
        slice_keys = [[split_val], ['Diagnostic'], antigen_list]
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
            if ci is not None:
                roc_suffix = '_'.join([suffix, 'ci'])
            else:
                roc_suffix = suffix
            #%%
            print('{} unique positive sera'.format(len(roc_df.loc[roc_df['serum type']=='positive', 'serum ID'].unique())))
            print('{} unique negative sera'.format(len(roc_df.loc[roc_df['serum type'] == 'negative', 'serum ID'].unique())))
            roc_plot_grid(roc_df, constants.RUN_PATH, '_'.join(['ROC', roc_suffix]), 'png', ci=ci, fpr=fpr, hue=hue)
    #%% Plot categorical scatter plot for episurvey
        if not cat_param_df.empty:
            sera_cat_list = cat_param_df['serum ID']
            slice_action = cat_param_df['serum ID action']
            hue = cat_param_df['hue']
            # plot specific slicing
            cat_df = slice_df(df_norm_sub, slice_action, 'serum ID', sera_cat_list)
            assert not cat_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
            sns.set_context("talk")
            g = sns.catplot(x="serum type", y="OD", hue=hue, col="antigen", kind="swarm",
                            palette=["r", "c", "y"], data=cat_df, col_wrap=3)
            plt.savefig(os.path.join(constants.RUN_PATH, 'catplot_{}.png'.format(suffix)),
                                          dpi=300, bbox_inches='tight')
            if fit_param_df['zoom']:
                g.set(ylim=(-0.05, 0.4))
                plt.savefig(os.path.join(constants.RUN_PATH, 'catplot_zoom_{}.png'.format(suffix)),
                                              dpi=300, bbox_inches='tight')
        #%% 4PL fit
        if not fit_param_df.empty:
            slice_action = fit_param_df['serum ID action']
            hue = fit_param_df['hue']
            dilution_df = slice_df(df_norm_sub, slice_action, 'serum ID', fit_param_df['serum ID'])
            dilution_df_fit = dilution_df.copy()
            dilution_df_fit = fit2df(dilution_df_fit, fourPL)
            sera_fit_list = dilution_df['serum ID'].unique()
            #%% plot standard curves
            sera_4pl_list = [' '.join([x, 'fit']) for x in sera_fit_list]
            antigens = dilution_df['antigen'].unique()
            markers = 'o'
            hue = 'serum ID'
            style = 'serum type'
            assert not dilution_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
            palette = sns.color_palette(n_colors=len(dilution_df[hue].unique()))
            print('plotting standard curves...')
            g = sns.lmplot(x="serum dilution", y="OD",
                            hue=hue, hue_order=sera_fit_list, col="antigen", ci='sd', palette=palette, markers=markers,
                             data=dilution_df, col_wrap=3, fit_reg=False, x_estimator=np.mean)
            palette = sns.color_palette(n_colors=len(dilution_df_fit[hue].unique()))
            for antigen, ax in zip(antigens, g.axes.flat):
                df_fit = dilution_df_fit[(dilution_df_fit['antigen'] == antigen)]
                sns.lineplot(x="serum dilution", y="OD", hue=hue, hue_order=sera_4pl_list, data=df_fit,
                             style=style, palette=palette,
                             ax=ax, legend=False)
                ax.set(xscale="log")
            plt.savefig(os.path.join(constants.RUN_PATH, 'fit_{}.jpg'.format(suffix)), dpi=300, bbox_inches='tight')
            if fit_param_df['zoom']:
                for antigen, ax in zip(antigens, g.axes.flat):
                    ax.set(ylim=[-0.05, 1.5])
                plt.savefig(os.path.join(constants.RUN_PATH, 'fit_{}_zoom.jpg'.format(suffix)), dpi=300, bbox_inches='tight')
        plt.close('all')
