import pandas as pd
import numpy as np
import os
import matplotlib
from matplotlib import pyplot as plt
from natsort import natsorted
import seaborn as sns
from interpretation.plotting import fourPL, fit2df, roc_plot_grid, joint_plot
from interpretation.report_reader import read_plate_info, read_antigen_info, read_pysero_output, read_scn_output, \
    slice_df, normalize_od, offset_od
import array_analyzer.extract.constants as constants

def read_config(input_dir):
    config_path = os.path.join(input_dir, constants.METADATA_FILE)
    if constants.METADATA_FILE not in os.listdir(input_dir):
        raise IOError("analysis config file not found, aborting")
    # check that the xlsx file contains necessary worksheets
    with pd.ExcelFile(config_path) as config_file:
        if not constants.LOAD_REPORT:
            assert ('pysero output dirs' in config_file.sheet_names) or \
            ('scienion output dirs' in config_file.sheet_names), \
            "sheet by name 'pysero output dirs' or 'scienion output dirs' are required " \
            "in analysis config file when load_report is False, aborting"

            if 'pysero output dirs' in config_file.sheet_names:
                ntl_dirs_df = pd.read_excel(config_file, sheet_name='pysero output dirs', index_col=0)
            if 'scienion output dirs' in config_file.sheet_names:
                scn_scn_df = pd.read_excel(config_file, sheet_name='scienion output dirs', index_col=0)
    return ntl_dirs_df, scn_scn_df


def analyze_od(input_dir, output_dir):
    sns.set_context("talk")
    font = {'size': 10, 'weight': 'normal', 'family': 'arial'}
    matplotlib.rc('font', **font)
    scn_scn_dirs = scn_scn_plate_ids = ntl_dirs = \
        ntl_slice_actions = ntl_well_ids = ntl_plate_ids = None
    ntl_dirs_df, scn_scn_df = read_config(input_dir)
    ntl_dirs = ntl_dirs_df['directory'].tolist()
    ntl_slice_actions = ntl_dirs_df['well action'].tolist()
    ntl_well_ids = ntl_dirs_df['well ID'].tolist()
    ntl_plate_ids = ntl_dirs_df['plate ID'].tolist()
    fig_path = os.path.join(output_dir, 'pysero_plots')

    sera_fit_list = ['Pool', 'mab', 'CR3022']
    sera_cat_list = ['Pool', 'mab', 'Blank', 'CR3022']
    sera_roc_list = sera_cat_list
    #%%
    os.makedirs(fig_path, exist_ok=True)
    load_master_report = True
    if not load_master_report:
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
        #%
        for data_folder, slice_action, well_id, plate_id in \
                zip(ntl_dirs, ntl_slice_actions, ntl_well_ids, ntl_plate_ids):
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

    stitched_pysero_df.loc[stitched_pysero_df['antigen']=='xIgG Fc', 'antigen type'] = 'Positive'
    #%% functions to compute ROC curves and AUC
    # for plate_id in stitched_pysero_df['plate_id'].unique():
    # for plate_id in ['plate_8']:
    #     slice_cols = ['pipeline', 'serum ID', 'plate_id']
    #     slice_keys = [['python'], sera_roc_list, [plate_id]]
    #     scn_psr_slice_actions = ['keep', 'drop', 'keep']

    norm_antigen = 'xIgG Fc'
    # norm_antigen = 'xkappa-biotin'
    # norm_antigen = None
    norm_group = 'plate'
    aggregate = 'mean'
    # aggregate = None
    # offset_antigen = 'GFP foldon'
    offset_antigen = None
    # norm_group = 'well'
    offset_group = 'well'
    antigen_list = ['SARS CoV2 N 50', 'SARS CoV2 RBD 250', 'SARS CoV2 spike 62.5']
    # sample_type = 'Orasure'
    sample_type = 'Serum'
    # slice_cols = ['serum ID', 'antigen type', 'antigen']
    # slice_keys = [sera_roc_list, ['Diagnostic'], antigen_list]
    # scn_psr_slice_actions = ['drop', 'keep', 'keep']
    slice_cols = ['serum ID', 'antigen type']
    slice_keys = [sera_roc_list, ['Diagnostic']]
    slice_actions = ['drop', 'keep', 'keep']
    fpr = 0.05
    # ci = 95
    ci = None
    hue = 'pipeline'
    for pipeline in stitched_pysero_df['pipeline'].unique():
    # for pipeline in ['nautilus']:
        df_norm = stitched_pysero_df.copy()
        df_norm = slice_df(df_norm, 'keep', 'pipeline', [pipeline])
        df_norm = slice_df(df_norm, 'keep', 'sample type', [sample_type])
        df_norm = normalize_od(df_norm, norm_antigen, group=norm_group)
        df_norm = offset_od(df_norm, offset_antigen, offset_group)
        suffix = '_'.join([pipeline, sample_type])
        if ci is not None:
            suffix = '_'.join([suffix, 'ci'])
        if norm_antigen is not None:
            suffix = '_'.join([suffix, norm_antigen, 'norm_per', norm_group])
        for col, action, key in zip(slice_cols, slice_actions, slice_keys):
            df_norm = slice_df(df_norm, action, col, key)
        roc_df = df_norm.copy()
        if aggregate is not None:
            roc_df = roc_df.groupby(['antigen', 'serum ID', 'well_id', 'plate_id', 'sample type',
                                     'serum type', 'serum dilution', 'pipeline', 'secondary ID',
                                     'secondary dilution'])['OD'].mean()
            roc_df = roc_df.reset_index()
            suffix = '_'.join([suffix, aggregate])
        # roc_df_split = roc_df.copy()
        # roc_df_split[['antigen', 'antigen conc']] = roc_df['antigen'].str.rsplit(n=1, expand=True)
        # roc_plot_grid(roc_df, fig_path, 'ROC')
        # roc_plot_grid(roc_df, fig_path, 'ROC_' + plate_id)
        roc_plot_grid(roc_df, fig_path, '_'.join(['ROC', suffix]), 'png', ci=ci, fpr=fpr, hue=hue)
        id_vars = [x for x in roc_df.columns if x not in ['False positive rate', 'True positive rate']]
        # roc_df = roc_df.melt(id_vars=id_vars,
        #                      var_name='category',
        #                      value_name='rate',
        #                      )
        # thr_plot_grid(roc_df, fig_path, '_'.join(['ROC_thr', suffix]), 'png')
    #%%
    df_serum = df_norm[['serum ID', 'serum type']].drop_duplicates()
    print(len(df_serum))
    print((df_serum['serum type']=='negative').sum())
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
    # norm_antigen = None
    norm_group = 'plate'
    # norm_group = 'well'
    # norm_antigen = 'xkappa-biotin'
    limit = [0, 2.0]
    neg_limit = [0, 0.1]
    x_cols = ['plate_3', 'plate_7', 'plate_4']
    y_cols = ['plate_9', 'plate_8', 'plate_10']
    # for pipeline in stitched_pysero_df['pipeline'].unique():
    for pipeline in ['nautilus']:
        df_norm = stitched_pysero_df.copy()
        df_norm = slice_df(df_norm, 'keep', 'pipeline', [pipeline])
        df_norm = normalize_od(df_norm, norm_antigen, group=norm_group)
        suffix = pipeline
        if norm_antigen is not None:
            suffix = '_'.join([pipeline, norm_antigen, 'norm_per', norm_group])

        pysero_df_pivot = pd.pivot_table(df_norm, values='OD',
                                     index=['well_id', 'antigen_row', 'antigen_col', 'serum ID', 'secondary ID', 'secondary dilution',
               'serum type', 'serum dilution', 'antigen', 'antigen type', 'pipeline'],
                                     columns=['plate_id'])
        pysero_df_pivot.reset_index(inplace=True)

        antigen_OD_df = slice_df(pysero_df_pivot, 'keep', 'antigen type', ['Diagnostic'])
        biotin_OD_df = slice_df(pysero_df_pivot, 'keep', 'antigen', ['xkappa-biotin'])
        # biotin_OD_df = slice_df(biotin_OD_df, 'keep', 'antigen type', ['Fiducial'])
        igg_OD_df = slice_df(pysero_df_pivot, 'keep', 'antigen', ['xIgG Fc'])
        antigen_pos_df = slice_df(antigen_OD_df, 'keep', 'serum type', ['positive'])
        antigen_neg_df = slice_df(antigen_OD_df, 'keep', 'serum type', ['negative'])

        for x_col, y_col in zip(x_cols, y_cols):
            # scatter_plot(pysero_df_pivot, x_col, y_col, fig_path, '_'.join(['OD_scatter', x_col, y_col, suffix]), xlim=limit, ylim=limit)
            # scatter_plot(antigen_OD_df, x_col, y_col, 'antigen', fig_path,
            #              '_'.join(['antigen_OD_scatter', x_col, y_col, suffix]), xlim=limit, ylim=limit)
            # scatter_plot(biotin_OD_df, x_col, y_col, 'biotin', fig_path,
            #              '_'.join(['biotin_OD_scatter', x_col, y_col, suffix]), xlim=limit, ylim=limit)
            # scatter_plot(igg_OD_df, x_col, y_col, 'igg', fig_path,
            #              '_'.join(['igg_OD_scatter', x_col, y_col, suffix]), xlim=limit, ylim=limit)
            joint_plot(antigen_OD_df, x_col, y_col, 'antigen type', 'antigen', fig_path,
                       '_'.join(['antigen_OD_joint', x_col, y_col, suffix]), bw='scott', n_levels=60, xlim=limit, ylim=limit)
            joint_plot(antigen_pos_df, x_col, y_col, 'antigen type', 'antigen', fig_path,
                       '_'.join(['antigen_pos_joint', x_col, y_col, suffix]), bw='scott', n_levels=60, xlim=limit, ylim=limit)
            joint_plot(antigen_neg_df, x_col, y_col, 'antigen type', 'antigen', fig_path,
                       '_'.join(['antigen_neg_joint', x_col, y_col, suffix]), bw='scott', n_levels=60, xlim=neg_limit, ylim=neg_limit)
        plt.close('all')

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
