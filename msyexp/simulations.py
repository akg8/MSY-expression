""" extracting and plotting simulation results """

import numpy as np 
import pandas as pd 
import scipy.stats as ss 
import matplotlib.pyplot as plt 
import seaborn as sns 

from . import xytools
from . import plotfuncs as pf

YORDER = ['AMELY', 'DDX3Y', 'EIF1AY', 'KDM5D', 'NLGN4Y', 'PRKY',
          'RPS4Y1', 'RPS4Y2', 'SRY', 'TBL1Y', 'TMSB4Y',
          'TXLNGY', 'USP9Y', 'UTY', 'ZFY',
          'PCDH11Y', 'TGIF2LY',
          'BPY2', 'CDY', 'DAZ', 'HSFY', 'PRY', 'RBMY', 'TSPY',
          'VCY','XKRY']

PAIRORDER = ['AMELY/AMELX', 'DDX3Y/DDX3X', 'EIF1AY/EIF1AX', 
             'KDM5D/KDM5C', 'NLGN4Y/NLGN4X', 'PRKY/PRKX',
             'RPS4Y1/RPS4X', 'SRY/SOX3', 'TBL1Y/TBL1X',
             'TMSB4Y/TMSB4X', 'USP9Y/USP9X',
             'UTY/KDM6A', 'ZFY/ZFX',
             'PCDH11Y/PCDH11X', 'TGIF2LY/TGIF2LX',
             'HSFY/HSFX', 'RBMY/RBMX', 'TSPY/TSPYL2',
             'VCY/VCX']

PALETTE_STRIP = {'cufflinks':'#225ea8', 'kallisto':'#d7301f', 
                 'fcts_unique':'#737373', 'salmon':'#238443'}
PALETTE_LIGHT = {'kallisto':'#fc8d59', 'cufflinks':'#41b6c4', 
                 'fcts_unique':'0.6', 'salmon':'#78c679'}
PALETTE_DARK = {'kallisto':'#990000', 'cufflinks':'#0c2c84', 
                'fcts_unique':'0.2', 'salmon':'#005a32'}

def extract_results(data, methods=['fcts_unique', 'cufflinks', 'kallisto'], 
                    log2ratio=False, yxratios=False):
    """ Return estatimed expression levels from relevant simulations 
    as log2(est/true) ratios or adjusted TPMs
    
    Parameters
    ----------
    data : DataFrame
        genes x expression level dataframe with 2-level hierarchical 
        column index: (method, simlib)
    methods : list
        methods to include in output
    log2ratio : bool
        extract results as log2(est / true); otherwise, extract as the
        adjusted TPM
    yxratios : bool
        return estimated Y/X ratios instead of expression levels;
        not compatible with log2ratio == True
    """
    res = {}
    for m in methods:
        if log2ratio:
            res[m] = np.log2((data[m] + 0.1) / (data['rsem_sim'] + 0.1))
        else:
            if not data['input'].loc[xytools.MSYFAMS].eq(0).all().all():
                scale = data['rsem_sim'] / data['input'] # adjustment for sampling in input
                res[m] = data[m] / scale
            else:
                res[m] = data[m]
    res = pd.concat(res, axis=1)

    if yxratios:
        dx = res.loc[xytools.XYPAIRS_X]
        dy = res.loc[xytools.XYPAIRS_Y]
        idx = ["{0}/{1}".format(gy, gx) for gx, gy in xytools.XYPAIRS]
        dx.index = dy.index = idx
        if log2ratio:
            res = dy - dx
        else:
            res = (dy + 0.1) / (dx + 0.1)

    return res

def summarize_results(data, methods=['fcts_unique', 'cufflinks', 'kallisto'],
                      yxratios=False):
    """ Return dataframe with summary statistics for simulation

    Parameters
    ----------
    data : DataFrame
        genes x expression level dataframe with 2-level hierarchical 
        column index: (method, simlib)
    methods : list
        methods to include in output
    """
    res = extract_results(data, log2ratio=False, methods=methods,
                          yxratios=yxratios)

    summ = {}
    
    if yxratios:
        _input_x = data.loc[xytools.XYPAIRS_X, 'input'].median(axis=1)
        _input_y = data.loc[xytools.XYPAIRS_Y, 'input'].median(axis=1)
        idx = ["{0}/{1}".format(gy, gx) for gx, gy in xytools.XYPAIRS]
        _input_x.index = _input_y.index = idx
        _input = pd.DataFrame({'all':(_input_y+0.1) / (_input_x+0.1)})

    else:
        res = res.loc[xytools.MSYFAMS]

        # input levels
        _input = data.loc[xytools.MSYFAMS, 'input'].median(axis=1)
        _input = pd.DataFrame({'all':_input})

    summ['input'] = _input

    for m in methods:
        res_m = res[m]

        # quantiles
        quant_m = res_m.quantile([0.025, 0.5, 0.975], axis=1).T
        quant_m.columns = ['2.5%', 'median', '97.5%']
        quant_m = quant_m.loc[:, ['median', '2.5%', '97.5%']]
        
        # unbiased
        quant_m['unbiased'] = (quant_m['2.5%'] < _input['all']) \
                                  & (quant_m['97.5%'] > _input['all'])
        
        # percent deviation
        abs_dev = np.abs(res_m.divide(_input['all'], axis=0)-1)
        quant_m['perc_dev(95%)'] = abs_dev.quantile(0.95, axis=1)
        summ[m] = quant_m

    summ = pd.concat(summ, axis=1)
    summ = summ.reindex(['input']+methods, level=0, axis=1)
    if yxratios:
        summ = summ.reindex(PAIRORDER, axis=0)
    else:
        summ = summ.reindex(YORDER, axis=0)
    return summ

def plot_results_boxplot(data, methods=['fcts_unique', 'cufflinks', 'kallisto'],
                         log2ratio=False, log2_thresh=-4, 
                         yxratios=False, log2yxratios=False,
                         palette_fc=PALETTE_LIGHT, palette_ec=PALETTE_DARK,
                         showfliers=False, showcaps=False):
    """ Boxplot of simulation results

    Parameters
    ----------
    data : DataFrame
        genes x expression level dataframe with 2-level hierarchical 
        column index: (method, simlib)
    methods : list
        methods to include in output
    log2ratio : bool
        plot results as log2(est / true)
    log2_thresh : float
        threshold for log2 values
    yxratios : bool
        plot Y/X ratios instead of Y levels
    log2yxratios : bool
        plot log2(estimated Y/X ratio)
    palette_fc : dict
        colors for box faces
    palette_ec : dict
        colors for box lines
    showfliers : bool
        show boxplot outliers
    showcaps : bool
        show boxplot caps
    """
    if (yxratios or log2yxratios):
        fig, ax = plt.subplots(figsize=(7.5, 3))
        fig.subplots_adjust(right=0.85, bottom=0.32, left=0.06, top=0.9)
    else:
        fig, ax = plt.subplots(figsize=(7.5, 2.5))
        fig.subplots_adjust(right=0.85, bottom=0.23, left=0.06, top=0.9)
    res0 = extract_results(data, log2ratio=log2ratio, methods=methods,
                           yxratios=(yxratios or log2yxratios))
    if yxratios:
        res0 = res0.loc[PAIRORDER]
    elif log2yxratios:
        res0 = np.log2(res0.loc[PAIRORDER])
    else:
        res0 = res0.loc[YORDER]
    
    for m in methods:
        res_m = res0.copy()
        methods0 = filter(lambda m0: m0 != m, methods)
        for m0 in methods0:
            res_m[m0] = np.nan
        fig, ax = _plot_sim_box(res_m, methods=methods, log2ratio=log2ratio,
                                log2_thresh=log2_thresh, yxratios=yxratios,
                                log2yxratios=log2yxratios,
                                palette=palette_fc, lcolor=palette_ec[m],
                                showfliers=showfliers, showcaps=showcaps,
                                ax=ax)        
    ax.spines['left'].set_linewidth(0.5)
    ax = pf.format_axis_labels(ax)
    return fig, ax

def _plot_sim_box(res, methods=['fcts_unique', 'cufflinks', 'kallisto'], 
                  log2ratio=False, log2_thresh=-4, yxratios=False, 
                  log2yxratios=False, showfliers=False, showcaps=False, 
                  ax=None, palette=PALETTE_LIGHT, lcolor='k'):
    """ auxiliary function for plotting a single hue of the boxplot
    """
    # format results
    if yxratios:
        xname = 'pair'
        res = res.reindex(PAIRORDER)
        if log2ratio:
            valname = 'Log2(est Y/X ratio / true Y/X ratio) [5th-95th percentile]'
            thresh = log2_thresh
        else:
            valname = 'Estimated Y TPM / X TPM [5th-95th percentile]'
            thresh = 0
            
    elif log2yxratios:
        xname = 'pair'
        res = res.reindex(PAIRORDER)
        valname = 'Log2(estimated Y TPM / X TPM) [5th-95th percentile]'
        thresh = log2_thresh
        
    else:
        xname = 'gene'
        res = res.reindex(YORDER)
        if log2ratio:
            valname = 'Log2(estimated TPM / true TPM) [5th-95th percentile]'
            thresh = log2_thresh
        else:
            valname = 'Estimated TPM [5th-95th percentile]'
            thresh = 0
    res.index.name = xname
    df = pd.melt(res.reset_index(), id_vars=xname, value_name=valname,
                 var_name=['method', 'libname'])
    df = df.loc[df.method.isin(methods)]
    df.loc[df[valname] < thresh, valname] = thresh
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 4))
    else:
        fig = ax.get_figure()
    ax = sns.boxplot(x=xname, y=valname, hue='method', data=df, whis=[5, 95],
                     hue_order=methods, linewidth=1,
                     palette=palette,
                     showcaps=showcaps, showfliers=showfliers,
                     boxprops={'edgecolor':lcolor, 'lw':0.5}, 
                     medianprops={'color':lcolor},
                     whiskerprops={'color':lcolor, 'lw':0.5}, 
                     capprops={'color':lcolor, 'lw':0.5},
                     flierprops={'markersize':2})
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
    for lab in ax.get_xticklabels():
        lab.set_rotation('vertical')
        lab.set_style('italic')
    for side in ['top', 'right', 'bottom']:
        ax.spines[side].set_visible(False)
    ax.grid(True, axis='y', ls=':')
    ax.set_xlabel('')
    return fig, ax

def plot_results_stripplot(data, methods=['fcts_unique', 'cufflinks', 'kallisto'],
                           log2ratio=False, log2_thresh=-4, 
                           yxratios=False, log2yxratios=False,
                           palette=PALETTE_STRIP):
    """ Stripplot of simulation results

    Parameters
    ----------
    data : DataFrame
        genes x expression level dataframe with 2-level hierarchical 
        column index: (method, simlib)
    methods : list
        methods to include in output
    log2ratio : bool
        plot results as log2(est / true)
    log2_thresh : float
        threshold for log2 values
    yxratios : bool
        plot Y/X ratios instead of Y levels
    log2yxratios : bool
        plot log2(estimated Y/X ratio)
    palette : dict
        colors for hue levels
    """
    if (yxratios or log2yxratios):
        fig, ax = plt.subplots(figsize=(7.5, 3))
        fig.subplots_adjust(right=0.85, bottom=0.32, left=0.06, top=0.9)
    else:
        fig, ax = plt.subplots(figsize=(7.5, 2.5))
        fig.subplots_adjust(right=0.85, bottom=0.23, left=0.06, top=0.9)
    res = extract_results(data, log2ratio=log2ratio, 
                          yxratios=(yxratios or log2yxratios))
    # format results
    if yxratios:
        xname = 'pair'
        res = res.reindex(PAIRORDER)
        if log2ratio:
            valname = 'log2(est Y/X ratio / true Y/X ratio)'
            thresh = log2_thresh
        else:
            valname = 'estimated Y TPM / X TPM'
            thresh = 0
            
    elif log2yxratios:
        xname = 'pair'
        res = np.log2(res.reindex(PAIRORDER))
        valname = 'log2(estimated Y TPM / X TPM)'
        thresh = log2_thresh
        
    else:
        xname = 'gene'
        res = res.reindex(YORDER)
        if log2ratio:
            valname = 'log2(estimated TPM / true TPM)'
            thresh = log2_thresh
        else:
            valname = 'estimated TPM'
            thresh = 0
    res.index.name = xname
    df = pd.melt(res.reset_index(), id_vars=xname, value_name=valname,
                 var_name=['method', 'libname'])
    df = df.loc[df.method.isin(methods)]
    df.loc[df[valname] < thresh, valname] = thresh

    ax = sns.stripplot(x=xname, y=valname, hue='method', data=df,
                       hue_order=methods, palette=palette,
                       alpha=0.4, dodge=True, s=2.5)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
    for lab in ax.get_xticklabels():
        lab.set_rotation('vertical')
        lab.set_style('italic')
    for side in ['top', 'right', 'bottom']:
        ax.spines[side].set_visible(False)
    ax.grid(True, axis='y', ls=':')
    ax.set_xlabel('')

    return fig, ax

def plot_sim_results(data, sim_name, log2yxratios=False, boxplot=True):
    """ Plot simulation results with default settings

    Parameters
    ----------
    data : DataFrame
        genes x expression level dataframe with 3-level hierarchical 
        column index: (sim_name, method, simlib)
    sim_name : str
        simulation name (top level of data.columns), e.g., 'sim_y1_x2'
    log2yxratios : bool
        plot log2(estimated Y/X); otherwise plot estimated Y TPM
    boxplot : bool
        plot boxplot; otherwise plot stripplot
    """
    if boxplot:
        fig, ax = plot_results_boxplot(data[sim_name], log2yxratios=log2yxratios)
    else:
        fig, ax = plot_results_stripplot(data[sim_name], 
                                         log2yxratios=log2yxratios)

    # simulation specific formatting
    if log2yxratios:
        if 'y5' in sim_name:
            title = 'Y/X ratios: Y = 5 TPM / X = 10 TPM'
        else:
            title = 'Y/X ratios: Y = 1 TPM / X = 2 TPM'
        ax.set_ylim([-4.2, 2.2])
        ax.spines['left'].set_bounds(-4, 2)
        ax.set_yticks([-4, -2, -1, 0, 2])
    else:
        if 'y5' in sim_name:
            yval = 5
            ax.set_ylim([-0.2, 10.2])
            ax.set_yticks([0, 4, 5, 6, 10])
            ax.spines['left'].set_bounds(0, 10)
        elif 'y1' in sim_name:
            yval = 1
            ax.set_ylim([-0.05, 2.05])
            ax.set_yticks(np.arange(0, 2.3, 0.5))
            ax.spines['left'].set_bounds(0, 2)
        else:
            yval = 0
            ax.set_ylim([-0.05, 1.05])
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
            ax.spines['left'].set_bounds(0, 1)
        title = 'Y genes ({} TPM)'.format(yval)
    ax.set_title(title)
    ax = pf.format_axis_labels(ax)
    ax = pf.format_spines(ax, bottom=False)

    return fig, ax

def plot_xy_uncorr(data, gx, gy, xy_connect=True, ax=None, 
                   showlegend=False, **kwargs):
    """ Plot results for X-Y uncorrelated simulations for one X-Y pair

    Note: code is hard-wired to assume an expected Y expression level of
    5 TPM with expected X expression levels in the range (0, 10)

    Parameters
    ----------
    data : DataFrame
        genes x expression level dataframe with 3-level hierarchical 
        column index: (sim_name, method, simlib)
    gx, gy : str
        X and Y homologs
    xy_connect : bool
        draw a connecting line between the X and Y estimated expression
        levels in each sample
    ax : axes_subplot
        existing axes to plot on
    showlegend : bool
        show the plot legend
    **kwargs : key, value pairs
        additional keyword arguments to pass to plt.scatter
    """
    ycolor = '#2166ac'
    ylcolor = '#4393c3'
    xcolor = '#d6604d'
    xlcolor = '#f4a582'
    
    res = extract_results(data['sim_y5_xuncorr'])
    x_in = data['sim_y5_xuncorr'].loc[gx, 'input']
    x_est = res.loc[gx, 'kallisto']
    y_est = res.loc[gy, 'kallisto']
    assert np.all(x_est.index == y_est.index)

    if ax is None:
        fig, ax = plt.subplots(figsize=(5, 5))
    else:
        fig = ax.get_figure()
    
    xlim = [-0.2, 10.2]
    ax.set_xlim(xlim)
    ax.plot([0, 10], [5, 5], ls='--', color=ycolor, lw=0.5, zorder=1)
    ax.plot([0, 10], [0, 10], ls='--', color=xcolor, lw=0.5, zorder=1)

    ax.scatter(x_in, x_est, edgecolors='none', facecolors=xcolor, 
               zorder=2, label=gx, **kwargs)
    ax.scatter(x_in, y_est, edgecolors='none', facecolors=ycolor, 
               zorder=2, label=gy, **kwargs)
    if xy_connect:
        for x, y0, y1 in zip(x_in, x_est, y_est):
            ax.plot([x, x], [y0, y1], color='0.6', ls='-', 
                    alpha=0.7, lw=0.5, zorder=1)
    
    if showlegend:
        ax.legend(loc='upper left', frameon=False)
    else:
        ax.set_title("{0}/{1}".format(gy, gx), style='italic', pad=0)
    ax.set_xlabel('Simulated X-homolog TPM', labelpad=2)
    ax.set_ylabel('Estimated TPM', labelpad=2)
    ax.set_xticks([0, 5, 10])
    ax.set_yticks([0, 5, 10])
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)
    for side in ['bottom', 'left']:
        ax.spines[side].set_color('k')
        ax.spines[side].set_linewidth(0.5)
        ax.spines[side].set_bounds(0, 10)

    # annotate with spearman correlation
    ok_idx = (~np.isnan(x_est)) & (~np.isnan(y_est))
    r, p = ss.spearmanr(x_est[ok_idx], y_est[ok_idx])
    sp_text = 'Spearman rho = {0:.3f}\np = {1:.3f}'.format(r, p)
    ax.text(0.2, 9.8, sp_text, size=6, color='k', 
            verticalalignment='top')
    
    return fig, ax

def plot_xy_uncorr_grid(data):
    """ Plot X-Y uncorrelated simulations for all X-Y pairs in a grid
    
    data : DataFrame
        genes x expression level dataframe with 3-level hierarchical 
        column index: (sim_name, method, simlib)
    """
    fig, axes = plt.subplots(5, 4, figsize=(7, 9))
    fig.subplots_adjust(left=0.07, right=0.97, top=0.96, bottom=0.07,
                        wspace=0.5, hspace=0.5)
    axes = axes.flatten()
    axes[-1].remove()
    xgenes = []
    ygenes = []
    for p in PAIRORDER:
        gy, gx = p.split('/')
        xgenes.append(gx)
        ygenes.append(gy)
    for gx, gy, ax in zip(xgenes, ygenes, axes[:19]):
        fig, ax = plot_xy_uncorr(data, gx, gy, ax=ax,
                                 s=8)
        ax = pf.format_axis_labels(ax)
    
    return fig, axes