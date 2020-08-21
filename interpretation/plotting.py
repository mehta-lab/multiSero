import itertools
import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from natsort import natsorted
from scipy import optimize as optimization
from sklearn.metrics import roc_curve, roc_auc_score


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
    s['False positive rate'], s['True positive rate'], s['threshold'] = \
        roc_curve(y_test, y_prob, pos_label='positive', drop_intermediate=False)
    try:
        s['AUC'] = [roc_auc_score(y_test, y_prob)] * len(s['False positive rate'])
    except ValueError as err:
        print('serum dilution {}, antigen {} only has {} serum type. {}'.
              format(df['serum dilution'].unique()[0], df['antigen'].unique()[0], y_test.unique()[0], err))
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
#%%

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
    ax.text(0.7 * xlim[1], 0.15 * ylim[1], 'Bias={:.3f}'.format(me), fontsize=16)
    ax.text(0.7 * xlim[1], 0.1 * ylim[1], 'Noise={:.3f}'.format(mae), fontsize=16)
    plt.savefig(os.path.join(output_path, ''.join([output_fname, '.jpg'])),
                dpi=300, bbox_inches='tight')
    plt.close()


def roc_plot(x, y, **kwargs):
    df = kwargs.pop('data')
    # y2 = kwargs.pop('y2')
    ax = plt.gca()
    df.plot(x=x, y=y, ax=ax, legend=False)
    # ax2 = ax.twinx()
    # color = 'g'
    # df.plot(x=x, y=y2, ax=ax2, color=color, style='--', legend=False)
    # ax2.set_ylabel(y2, color=color)  # we already handled the x-label with ax1
    # ax2.tick_params(axis='y', labelcolor=color)


def roc_plot_grid(roc_df, fig_path, fig_name, ext, col_wrap=3):
    # Plot ROC curves
    hue = 'serum dilution'
    fpr = 0.01
    antigens = natsorted(roc_df['antigen'].unique())
    sns.set_context("notebook")
    assert not roc_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
    palette = sns.color_palette(n_colors=len(roc_df[hue].unique()))
    print('plotting ROC curves...')
    g = sns.FacetGrid(roc_df, hue=hue, col="antigen", col_order=antigens, col_wrap=col_wrap, aspect=1,
                      xlim=(-0.05, 1), ylim=(0, 1.05))
                      # hue_kws={'linestyle': ['-', '--', '-.', ':']})
    g = (g.map_dataframe(roc_plot, 'False positive rate', 'True positive rate', y2='threshold'))
    for antigen, ax in zip(antigens, g.axes.flat):
        sub_df = roc_df[roc_df['antigen'] == antigen]
        tpr = np.interp(fpr, sub_df['False positive rate'], sub_df['True positive rate'])
        ax.plot([fpr, fpr], [0, tpr], linewidth=1, color='g', linestyle='--', alpha=1)
        ax.plot([-0.05, fpr], [tpr, tpr], linewidth=1, color='g', linestyle='--', alpha=1)
        auc = sub_df['AUC'].unique()[0]
        ax.set_title(antigen)
        ax.text(0.6, 0.15, 'AUC={:.3f}'.format(auc), fontsize=12)  # add text
        ax.text(fpr + 0.05, tpr - 0.2, 'sensitivity={:.3f}\nspecificity={:.3f}'.format(tpr, 1-fpr),
                fontsize=12, color='g')  # add text
    plt.savefig(os.path.join(fig_path, '.'.join([fig_name, ext])),
                             dpi=300, bbox_inches='tight')
    plt.close()


def thr_plot_grid(roc_df, fig_path, fig_name, ext, col_wrap=3):
    # Plot ROC curves
    hue = 'category'
    antigens = natsorted(roc_df['antigen'].unique())
    sns.set_context("notebook")
    assert not roc_df.empty, 'Plotting dataframe is empty. Please check the plotting keys'
    palette = sns.color_palette(n_colors=len(roc_df[hue].unique()))
    print('plotting ROC curves...')
    g = sns.FacetGrid(roc_df, hue=hue, col="antigen", col_order=antigens, col_wrap=col_wrap, aspect=1,
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
    plt.savefig(os.path.join(fig_path, '.'.join([fig_name, ext])),
                             dpi=300, bbox_inches='tight')
    plt.close()