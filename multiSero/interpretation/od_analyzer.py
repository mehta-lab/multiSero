import pandas as pd
import os
import re
import itertools
from matplotlib import pyplot as plt
import seaborn as sns
from matplotlib.transforms import Bbox

from multiSero.interpretation.plotting import roc_plot_grid, total_plots
from multiSero.interpretation.report_reader import slice_df, normalize_od, read_output_batch
import multiSero.array_analyzer.extract.constants as constants


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
    antigen_list = []
    antigen_list.append(plot_setting_df['antigens to plot']) #define list with first element, then the rest, can be one line
    antigen_list.append('SARS CoV2 RBD 250')
    antigen_list.append('SARS CoV2 spike 62.5')

    suffix = ''
    df_norm = normalize_od(stitched_multisero_df.copy(), norm_antigen, group=norm_group)
    if norm_antigen is not None:
        suffix = '_'.join([suffix, norm_antigen, 'norm_per', norm_group])

    if aggregate is not None:
        df_norm = df_norm.groupby(['antigen', 'antigen type', 'serum ID', 'well_id', 'plate ID', 'sample type',
                                 'serum type', 'serum dilution', 'serum cat', 'pipeline', 'secondary ID',
                                 'secondary dilution', 'visit value'])['OD'].mean().reset_index()
        suffix = '_'.join([suffix, aggregate])

    for split_val in split_plots_vals:
        split_suffix = suffix
        if split_val is not None:
            split_suffix = '_'.join([split_suffix, split_val])
        df_norm_sub = slice_df(df_norm, 'keep', split_plots_by, split_val)
        slice_cols = [split_plots_by, 'antigen type', 'antigen']
        slice_keys = [[split_val], ['Diagnostic'], antigen_list]
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
            #multisero_df #make sure that multisero_df is tidy
            ms_tidy_df = df_norm_sub[['antigen','OD','visit value','serum cat','serum type','serum ID','serum dilution']]
            ms_tidy_df.drop(ms_tidy_df.loc[ms_tidy_df['serum dilution'] > 0.00005].index, inplace=True)
            ms_tidy_df.drop(ms_tidy_df.loc[(ms_tidy_df['serum cat'] != '[\'COVID-Vax+\']')
                                           & (ms_tidy_df['serum cat'] != '[\'COVID+Vax+\']')
                                           & (ms_tidy_df['serum cat'] != '[\'COVID+Vax-\']')].index, inplace=True)
            #ms_tidy_df.drop(ms_tidy_df.loc[(ms_tidy_df['serum cat'] != 'COVID+Vax+') & (ms_tidy_df['serum cat'] != 'COVID+Vax-')
                                           #& (ms_tidy_df['serum cat'] != 'COVID-Vax+')].index, inplace=True)
            #for each in ms_tidy_df['serum type'].unique():

            deltaod = ms_tidy_df.copy(deep=True)
            deltaod.set_index(['serum type','antigen'],inplace=True)
            #grouped = deltaod.groupby(['serum type','antigen'])['OD'].mean()#still need to separate by visit val
            #grouped = deltaod.groupby(['serum type', 'antigen','visit value'])['OD'] #not complete
            #abc = deltaod.loc[("serum type","antigen"),"visit value"]
            #deltaod.groupby(['serum type', 'antigen','visit value'])['OD'].transform(lambda x: x[0] / x[1])
            #deltaod.groupby(level=0)[3].transform(lambda x: x[0] / x[1])

            #for antigen, new_df in deltaod.groupby(level=0):
                #print(new_df)
            group = ms_tidy_df.groupby('serum type')
            #df3 = group.apply(lambda x: x['OD'])
            df2 = group.apply(lambda x: x['visit value'].unique())
            #find mean OD for patients for each serum ID where time bin = early
            for element in df2.index:
                #if element == ms_tidy_df['serum type']
                if (len(df2[element]) > 1):
                    #ms_tidy_df['dOD'] = deltaod.groupby(['serum type', 'antigen', 'visit value'])['OD'].transform(lambda x: x[0] / x[1])
                    if (df2[element][0][2:-2] > df2[element][1][2:-2]):
                        ms_tidy_df.loc[ms_tidy_df['visit value'] == df2[element][0], 'vaccine availability'] = 'post-vax'
                        ms_tidy_df.loc[ms_tidy_df['visit value'] == df2[element][1], 'vaccine availability'] = 'pre-vax'
                    if (df2[element][0][2:-2] < df2[element][1][2:-2]):
                        ms_tidy_df.loc[ms_tidy_df['visit value'] == df2[element][0], 'vaccine availability'] = 'pre-vax'
                        ms_tidy_df.loc[ms_tidy_df['visit value'] == df2[element][1], 'vaccine availability'] = 'post-vax'
                else:
                    print(element)
            #if ms_tidy_df['visit value']

            prevaxdf = ms_tidy_df.loc[ms_tidy_df['vaccine availability'] == 'pre-vax']
            prevaxdf['old OD'] = prevaxdf['OD']
            #prevaxdf.set_index(['serum type', 'antigen'], inplace=True)
            postvaxdf = ms_tidy_df.loc[ms_tidy_df['vaccine availability'] == 'post-vax']
            postvaxdf['new OD'] = postvaxdf['OD']
            #prevaxdf.set_index(['serum type', 'antigen'], inplace=True)
            result = pd.merge(prevaxdf, postvaxdf, how='inner', on=['antigen', 'serum type'])
            result['delta OD'] = result['new OD'] - result['old OD']
            """
            #delta_od_df =
            #THE FOLLOWING DOESN'T MAKE SENSE BECAUSE YOU'RE NOT SELECTING FOR ANTIGEN TYPES.
            for each in ms_tidy_df['serum type'].unique(): #TO MAKE IT MAKE SENSE I THINK YOU NEED TO DO IT LIKE THEY DO IT IN PLOTTING.PY
                ms_tidy_df.loc[((ms_tidy_df['serum type'] == each) & (ms_tidy_df['time bin'] == 'early')),
                               'mean early OD'] = ms_tidy_df['OD'].loc[(ms_tidy_df['serum type'] == each) &
                                                                       (ms_tidy_df['time bin'] == 'early')].mean()
                ms_tidy_df.loc[((ms_tidy_df['serum type'] == each) & (ms_tidy_df['time bin'] == 'late')),
                               'mean late OD'] = ms_tidy_df['OD'].loc[(ms_tidy_df['serum type'] == each) &
                                                                       (ms_tidy_df['time bin'] == 'late')].mean()
            """

                #delta_od_df['serum type'] = each
                #delta_od_df['delta OD'] = ms_tidy_df['mean late OD'].loc[ms_tidy_df['serum type'] == each] - ms_tidy_df['mean early OD'].loc[ms_tidy_df['serum type'] == each]
            #
            #
            sera_cat_list = cat_param_df['serum ID']
            slice_action = cat_param_df['serum ID action']
            split_subplots_by = cat_param_df['split subplots by']
            hue = cat_param_df['hue']
            # plot specific slicing
            #cat_df = slice_df(df_norm_sub, slice_action, 'serum ID', sera_cat_list) #serum ID --> antigen
            ms_tidy_df.loc[ms_tidy_df['serum cat'] == '[\'COVID-Vax+\']', 'vaccine availability'] = 'post-vax'
            #ms_tidy_df.drop(ms_tidy_df.loc[ms_tidy_df['serum cat'] == 'neg'].index, inplace=True)
            cat_df = slice_df(ms_tidy_df, slice_action, 'serum ID', sera_cat_list)
            assert not cat_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
            sns.set_context("talk")

            ## FIGURE 5A -- CATPLOT OF OD -- VIOLIN NOT INCLUDED
            fig_palette = ["#ee2cef", "#21e020", "#248ff9", "#e7e5e6"]
            g = sns.catplot(x="vaccine availability", y="OD", hue=hue, col=split_subplots_by, kind="swarm", legend=False,
                            palette=fig_palette, order=["pre-vax", "post-vax"], dodge=True, data=cat_df, col_wrap=3)
            #g.map_dataframe(sns.violinplot, x="vaccine availability", y="OD", hue=hue, color="0.8", dodge=True,
                            #order=["pre-vax", "post-vax"], alpha=0.3)
            #g.map_dataframe(sns.boxplot, x="vaccine availability", y="OD", medianprops={'color': 'k', 'ls': '-', 'lw': 2},
                            #whiskerprops={'visible': False}, showfliers=False, showbox=False, showcaps=False, zorder=10,
                            #dodge=True)
            #plt.legend(bbox_to_anchor=(1.02, 0.5), loc='center left', borderaxespad=0)
            #TODO: add median OD value labels per antigen per cat group
            #DONE, opted for boxplot instead
            from numpy import median
            median_palette = ["#a318a4","#185ea4","#18a418"]
            #clr = itertools.cycle(median_palette

            g.map_dataframe(sns.boxplot, x="vaccine availability", y="OD", hue=hue,
                            #medianprops={'color': 'k', 'ls': '-', 'lw': 3},
                            whiskerprops={'visible': False},
                            showfliers=False, showbox=False, showcaps=False, zorder=10, order=["pre-vax", "post-vax"],
                            palette=median_palette,
                            dodge=0.55) #change hue

            for ax in g.axes.flat:
                #print(range(15)[1::2])
                for line in ax.get_lines()[2:10:2]:
                    line.set_color(median_palette[0])
                for line in ax.get_lines()[10::2]:
                    line.set_color(median_palette[1])
                for line in ax.get_lines()[1::2]:
                    line.set_color(median_palette[2])
            #g.map_dataframe(sns.pointplot, x="vaccine availability", y="OD", hue=hue, dodge=0.55, join=False, ci=None,
                            #estimator=median, order=["pre-vax", "post-vax"], zorder=10, markers="*",palette=fig_palette)
            """
            hue_labels = ['COVID+/Vax+', 'COVID+/Vax-','COVID-/Vax+','median']
            g.add_legend(legend_data={
                key: value for key, value in zip(hue_labels, g._legend_data.values())})
            """
            ###
            #legend_labels, _ = g._legend_data.keys()
            #ax.legend(legend_labels, ['man1', 'woman1', 'child1'], bbox_to_anchor=(1, 1))
            ###
            q = plt.legend()
            q = plt.legend(bbox_to_anchor=(1.02, 0.5), loc='center left', borderaxespad=0, handleheight=0.05,
                           edgecolor="#000000")
            q.get_texts()[0].set_text('COVID+/Vax+ median (N:;RBD:;Spike:)')
            q.get_texts()[1].set_text('COVID+/Vax- median (N:;RBD:;Spike:)')
            q.get_texts()[2].set_text('COVID+/Vax+ median (N:;RBD:;Spike:)')
            #legls = q.legendHandles[0:3]
            """
            handlelength = 1, handleheight = 1
            
            for lh in legls:
                lh.set_alpha(0.1) #control tickness through bbox_inches attribute
                print(lh.figure.bbox_inches)
                lh.figure.bbox_inches = Bbox(x0=0.0,y0=0.0,x1=15.0,y1=2.0,points=[])
                #lh._sizes = [2]
            """


            #plt.xlabel("vaccine availability")
            #plt.ylabel("OD")
            g.set_axis_labels("vaccine availability", "OD")
            g.set_xticklabels(rotation=0, horizontalalignment='center')

            plt.savefig(os.path.join(constants.RUN_PATH, 'catplot_{}.png'.format(split_suffix)),
                        dpi=300, bbox_inches='tight')
            if cat_param_df['zoom']:
                g.set(ylim=(-0.05, 0.4))
                plt.savefig(os.path.join(constants.RUN_PATH, 'catplot_zoom_{}.png'.format(split_suffix)),
                            dpi=300, bbox_inches='tight')

            ## FIGURE 5B -- PLOTTING DELTA OD -- VIOLIN NOT INCLUDED -- MEDIAN LINE ADDED
            hue = 'serum cat_x'
            g = sns.catplot(x="serum cat_x", y="delta OD", hue=hue, col=split_subplots_by, kind="swarm", palette= fig_palette,
                            data=result, col_wrap=3, legend=False)
            #g.map_dataframe(sns.violinplot, x="serum cat_x", y="delta OD", color="0.8", hue=hue, alpha=0.3, dodge=False)
            g.map_dataframe(sns.boxplot, x="serum cat_x", y="delta OD", medianprops={'color': 'k', 'ls': '-', 'lw': 2},
                            whiskerprops={'visible': False}, showfliers=False, showbox=False, showcaps=False, zorder=10,
                            dodge=False) #want to make the width a little less wide and maybe change line style
            #medians = result.groupby(['serum cat_x'])['delta OD'].mean()

            #plt.legend(bbox_to_anchor=(1.02, 0.5), loc='center left', borderaxespad=0)
            hue_labels = ['COVID+/Vax+', 'COVID+/Vax-','median']
            g.add_legend(legend_data={
                key: value for key, value in zip(hue_labels, g._legend_data.values())})
            g.set_axis_labels("cohort", "ΔOD")
            ls = ['COVID+/Vax+','COVID+/Vax-']
            g.set_xticklabels(ls, rotation=0)
            #plt.ylabel("ΔOD")
            plt.savefig(os.path.join(constants.RUN_PATH, 'catplot_deltaod_{}.png'.format(split_suffix)),
                        dpi=300, bbox_inches='tight')
            if cat_param_df['zoom']:
                g.set(ylim=(-0.05, 0.4))
                plt.savefig(os.path.join(constants.RUN_PATH, 'catplot_od_zoom_{}.png'.format(split_suffix)),
                            dpi=300, bbox_inches='tight')
        #%% 4PL fit
        if not fit_param_df.empty:
            slice_action = fit_param_df['serum ID action']
            hue = fit_param_df['hue']
            dilution_df = slice_df(df_norm_sub, slice_action, 'serum ID', fit_param_df['serum ID'])
            split_subplots_by = fit_param_df['split subplots by']
            #total_plots(dilution_df, constants.RUN_PATH, 'fit_{}'.format(split_suffix), 'png', hue=hue,
                        #zoom=fit_param_df['zoom'], split_subplots_by=split_subplots_by, col_wrap=2)
        plt.close('all')

