""" functions for analyzing Y/X epxression ratios """
import numpy as np 
import pandas as pd 
import scipy.stats as ss 
import statsmodels.api as sm 
import matplotlib.pyplot as plt
import seaborn as sns 

from . import xytools
from . import general
from . import plotfuncs as pf 

def calc_ratios(data, genes1=xytools.XYPAIRS_Y, genes2=xytools.XYPAIRS_X, 
			    min_exp=1., pseudo=0.5):
    """ Calculate sample-level log2-expression ratios of genes

    Parameters
    ----------
    data : DataFrame
        genes x observations dataframe
    genes1 : list-like
        first set of genes (numerators); or if genes2 is None, 
        this list contains length-2 list-like objects where the
        first entry in each is the numerator gene and the second
        is the denominator gene, e.g. [('DDX3X', 'DDX3Y'), ...]
    genes2 : list-like
        second set of genes (denominators)
    min_exp : float
        only calculate a ratio in observations where at least
        one of the two genes has an expression level at least
        this great
    pseudo : float
        pseudocounts to add to numerator and denominator prior
        to taking log2 ratio
    """
    if genes2 is None:
        _genes1 = []
        _genes2 = []
        for g1, g2 in genes1:
            _genes1.append(g1)
            _genes2.append(g2)
        genes1 = _genes1
        genes2 = _genes2
    else:
        if len(genes1) != len(genes2):
            raise ValueError("gene lists must be the same length")
    if not set(genes1).issubset(data.index):
        raise ValueError("some genes not found in data.index")
    if not set(genes2).issubset(data.index):
        raise ValueError("some genes not found in data.index")

    d1 = data.loc[genes1]
    d2 = data.loc[genes2]
    idx = ["{0}/{1}".format(g1, g2) for g1, g2 in zip(genes1, genes2)]
    d1.index = d2.index = idx
    ratio = np.log2((d1 + pseudo) / (d2 + pseudo))
    ratio[(d1 < min_exp) & (d2 < min_exp)] = np.nan 
    return ratio

def calc_ratios_by_tissue(data, meta, genes1=xytools.XYPAIRS_Y, 
                          genes2=xytools.XYPAIRS_X, min_samp_exp=1.0, 
                          min_tissue_exp=1.0, pseudo=0.5, func='median', 
                          tissuecol='TISSUE'):
    """ Calculate median expression ratios in each tissue

    Parameters
    ----------
    data : DataFrame
        genes x samples dataframe of TPMs for individual genes
    meta : DataFrame
        metadata dataframe
    genes1 : list-like
        first set of genes (numerators); or if genes2 is None, 
        this list contains length-2 list-like objects where the
        first entry in each is the numerator gene and the second
        is the denominator gene, e.g. [('DDX3X', 'DDX3Y'), ...]
    genes2 : list-like
        second set of genes (denominators)
    min_samp_exp : float
        only calculate a ratio in samples where at least
        one of the two genes has an expression level at least
        this great
    min_tissue_exp : float
        only report a pair's expression ratio in a given tissue
        if at least one of the genes has a median expression level
        greater than or equal to this value
    pseudo : float
        pseudocounts to add to numerator and denominator prior
        to taking log2 ratio
    func : {'median'|'mean'}
        function for aggregating ratios in each tissue
    """
    if genes2 is None:
        _genes1 = []
        _genes2 = []
        for g1, g2 in genes1:
            _genes1.append(g1)
            _genes2.append(g2)
        genes1 = _genes1
        genes2 = _genes2
    ratios = calc_ratios(data, genes1, genes2, min_exp=min_samp_exp, 
                         pseudo=pseudo)
    rbt = general.stat_by_tissue(ratios, meta, func=func, tissuecol=tissuecol)
    
    # get tissue expression levels
    mbt = general.stat_by_tissue(data.loc[genes1+genes2], meta, 
                                 func=func, tissuecol=tissuecol)
    mbt1 = mbt.loc[genes1]
    mbt2 = mbt.loc[genes2]
    mbt1.index = mbt2.index = ratios.index
    rbt[(mbt1 < min_tissue_exp) & (mbt2 < min_tissue_exp)] = np.nan
    return rbt

def yxratio_sig_test(data, meta, min_samps=10, min_exp=1., verbose=False,
                     alpha=0.01, genes_x=xytools.XYPAIRS_X,
                     genes_y=xytools.XYPAIRS_Y):
    """ Perform Wilcoxon rank-sum tests for Y/X ratio != 0
    
    Parameters
    ----------
    data : pd.DataFrame
        genes x samples dataframe; X and Y homologs will be extracted
    meta : pd.DataFrame
        samples metadata
    min_samps : int
        require that the X or Y homolog of an X-Y pair is expressed 
        in at least this many samples from a tissue to perform a 
        significance test
    min_exp : float
        threshold for expression
    alpha : float
        critical value (as applied to q-values) for determining whether
        X=Y, X>Y, or Y>X.

    Returns
    -------
    res : DataFrame
        longform data frame with log2(Y/X ratios), p-values, and 
        BH-adjusted p-values (i.e., q-values) for each gene-pair/tissue
        combination
    """
    pseudo = 0.5
    meta_ = meta.set_index('SAMPID')
    rbt = calc_ratios_by_tissue(data, meta, min_samp_exp=min_exp,
                                min_tissue_exp=min_exp, pseudo=pseudo,
                                genes1=genes_y, genes2=genes_x)
    
    d_x = data.loc[genes_x].copy()
    d_y = data.loc[genes_y].copy()
    
    pvals = pd.DataFrame(np.empty_like(rbt), index=rbt.index, 
                         columns=rbt.columns)

    for i, t in enumerate(rbt.columns):
        if verbose:
            general._print("{0} [{1}/{2}]...".format(t, i+1, rbt.shape[1]))
        # extract tissue data
        dt = data.loc[:, meta_.loc[data.columns, 'TISSUE']==t].copy()
        d_xt = dt.loc[genes_x].copy()
        d_yt = dt.loc[genes_y].copy()
        idx = ["{0}/{1}".format(y, x) for y, x in zip(genes_y, genes_x)]
        d_xt.index = d_yt.index = idx

        for p in rbt.index:
            if np.isnan(rbt.at[p, t]):
                pvals.at[p, t] = np.nan
            else:
                d_xt_p = d_xt.loc[p].copy()
                d_yt_p = d_yt.loc[p].copy()
                # check to make sure there are enough samples for test
                oksamps = (d_xt_p > min_exp) | (d_yt_p > min_exp)
                if oksamps.sum() < min_samps:
                    pvals.at[p, t] = np.nan
                else:
                    d_xt_p = d_xt_p.loc[oksamps]
                    d_yt_p = d_yt_p.loc[oksamps]

                    _, pvals.at[p, t] = ss.wilcoxon(d_yt_p, d_xt_p)
    
    # combine with log2(Y/X ratios)
    pvals.index.name = 'pair'
    pvals = pd.melt(pvals.reset_index(), id_vars='pair', 
                    var_name='tissue', value_name='pval')
    rbt.index.name = 'pair'
    rbt = pd.melt(rbt.reset_index(), id_vars='pair', 
                  var_name='tissue', value_name='log2(Y/X ratio)')
    rbt['Y/X ratio'] = 2**rbt['log2(Y/X ratio)']
    res = pd.merge(rbt, pvals, on=['pair', 'tissue'])
    
    # FDR
    _df = res[['pair', 'tissue', 'pval']].copy()
    _df = _df.loc[_df.notnull().all(axis=1)]
    _, _df['pval_adj'], _, _ = sm.stats.multipletests(_df['pval'], 
                                                      method='fdr_bh')
    _ = _df.pop('pval')
    res = pd.merge(res, _df, on=['pair', 'tissue'], how='outer')

    # annotate result
    res['status'] = 'X = Y'
    res.loc[(res.pval_adj < alpha) & (res['Y/X ratio'] > 1), 'status'] = 'Y > X'
    res.loc[(res.pval_adj < alpha) & (res['Y/X ratio'] < 1), 'status'] = 'X > Y'
    res.loc[res.pval.isnull(), 'status'] = np.nan
    
    return res

def plot_yxratios(res, color_by_tissue=False, statcol='status', 
                  min_l2ratio=-4, sort_genes=True):
    """ plot stripplot with Y/X expression ratios

    res : DataFrame
        as output by yxratio_sig_test
    color_by_tissue : bool
        color points by testis/non-testis; otherwise color by significance
    statcol : 'fdr_bh' or 'bonferroni'
        the column of res to use when coloring by Y>X/X>Y significance
    min_l2ratio : float
        trim low log2 Y/X ratios at this value 
    sort_genes : bool
        sort genes by median Y/X ratio
    """
    df = res.copy()
    palette = {'Y > X':'#045390', 'X > Y':'#da7b28', 'X = Y':'#a7a9ac',
           'testis':'r', 'other':'0.2'}

    if sort_genes:
        p_order = df.groupby('pair').median()['Y/X ratio']
        p_order = list(p_order.sort_values(ascending=False).index)

    df['tissue_'] = 'other'
    df.loc[df.tissue=='Testis', 'tissue_'] = 'testis'
    df.loc[df['log2(Y/X ratio)'] < min_l2ratio, 'log2(Y/X ratio)'] = min_l2ratio
    if color_by_tissue:
        hue = 'tissue_'
    else:
        hue = statcol

    fig, ax = plt.subplots(figsize=(4, 2.5))
    fig.subplots_adjust(bottom=0.347, right=0.694, top=0.920, left=0.108)
    ax = sns.stripplot(x='pair', y='log2(Y/X ratio)', hue=hue, data=df,
                       order=p_order, alpha=0.5, palette=palette, s=3)

    # plot highlighted point for RPS4X/Y in skeletal muscle
    x_ = p_order.index('RPS4Y1/RPS4X')
    y_ = res.loc[(res['pair']=='RPS4Y1/RPS4X') & (res['tissue']=='Skeletal Muscle'), 'log2(Y/X ratio)'].iloc[0]
    ax.scatter([x_], [y_], facecolors='#e92e2f', edgecolors='y', 
               zorder=3, s=10)

    for lab in ax.get_xticklabels():
        lab.set_rotation('vertical')
        lab.set_style('italic')
    ax.grid(True, axis='y', ls=':', lw=0.5, color='0.7')
    yticks = np.array([-4, -2, -1, 0, 1, 2])
    ytlabs = list(map(lambda y: "{:.2f}".format(y), 2.**yticks))
    ax.set_yticks(yticks)
    ax.set_ylim([-4, 3])
    ax.set_yticklabels(ytlabs)
    xlim = [-0.5, len(p_order)-0.5]
    ax.set_xlim(xlim)
    ax.plot(xlim, [0, 0], lw=0.5, ls='-', color='0.5')
    
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False,
              title='status\n(FDR < 1%)', title_fontsize=7,
              fontsize=7, markerscale=0.5)
    ax.set_ylabel('Y/X expression ratio in tissue')
    for lab in ax.get_xticklabels():
        lab.set_style('italic')
    
    ax.spines['left'].set_bounds(-4, 3)
    ax = pf.format_spines(ax, bottom=False)
    ax = pf.format_axis_labels(ax, labsize=7, ticklabsize=6)
    ax.set_xlabel('')
    return fig, ax

def compare_gtex_hpa_yxratios(gratio, hratio, gmeta, hmeta, tissue,
                              broad_only=False, show_hpa_points=False,
                              show_gtex_ci=False, show_hpa_ci=False,
                              smallfig=False,
                              calc_mwu_pvals=False, p_thresh=0.05,
                              annot=False):
    """ Plot Y/X ratios from GTEx and HPA for one tissue 
    
    Parameters
    ----------
    gratio, hratio : DataFrame
        pair x log2(Y/X ratio) dataframes from GTEx, HPA
    gmeta, hmeta : DataFrame
        metadata DataFrames from GTEx, HPA
    tissue : str
        tissue to compare
    broad_only : bool
        use only broadly expressed genes
    """
    if broad_only:
        pairs = ['{0}/{1}'.format(gy, gx) for gy, gx in zip(xytools.XYPAIRS_Y9, xytools.XYPAIRS_X9)]
        gratio = gratio.loc[gratio.index.isin(pairs)]
        hratio = hratio.loc[hratio.index.isin(pairs)]
    else:
        pairs = [p for p in gratio.index]
    
    # extract tissue-level estimates
    hmeta_ = hmeta.set_index('SAMPID')
    gmeta_ = gmeta.set_index('SAMPID')
    hvals = hratio.loc[:, hmeta_.loc[hratio.columns, 'TISSUE']==tissue]
    gvals = gratio.loc[:, gmeta_.loc[gratio.columns, 'TISSUE']==tissue]
    _df = pd.DataFrame({'gtex':gvals.mean(axis=1), 
                        'hpa':hvals.mean(axis=1)})
    _df = _df.loc[_df.notnull().all(axis=1)]
    
    if smallfig:
        marker_s = 20
        lw = 0.5
    else:
        marker_s = 40
        lw = 1
    fig, ax = plt.subplots(figsize=(4, 4))
    fig.subplots_adjust(left=0.2, bottom=0.2, right=0.95)
    if calc_mwu_pvals:
        mwu_pvals = {}
        for p in pairs:
            u, pval = general.mannwhitney_permute(hvals.loc[p].dropna(),
                                                  gvals.loc[p].dropna())
            mwu_pvals[p] = pval
        _df['pval'] = pd.Series(mwu_pvals)
        ax.scatter(_df.loc[_df['pval'] >= p_thresh, 'gtex'],
                   _df.loc[_df['pval'] >= p_thresh, 'hpa'],
                   edgecolors='none', facecolors='k', alpha=0.7,
                   s=marker_s, zorder=2, label='P >= {:.2e}'.format(p_thresh))
        ax.scatter(_df.loc[_df['pval'] < p_thresh, 'gtex'],
                   _df.loc[_df['pval'] < p_thresh, 'hpa'],
                   edgecolors='none', facecolors='r', alpha=0.7,
                   s=marker_s, zorder=2, label='P < {:.2e}'.format(p_thresh))
        if (_df['pval'] < p_thresh).any():
            for pair, row in _df.loc[_df['pval'] < p_thresh].iterrows():
                ax.text(row['gtex'], row['hpa'], pair, size=6,
                        color='r', horizontalalignment='left',
                        verticalalignment='top')
        ax.legend(title='MW-U permute')

    else:
        ax.scatter(_df.gtex, _df.hpa, edgecolors='none',
                   facecolors='k', alpha=0.7, s=marker_s, zorder=2)
        if annot:
            for pair, row in _df.iterrows():
                ax.text(row['gtex'], row['hpa'], pair, size=6,
                        color='r', horizontalalignment='left',
                        verticalalignment='top')
    
    # count number of samples per tissue for plot labels
    gcounts = gvals.shape[1]  # gmeta_.loc[gratio.columns, 'TISSUE'].value_counts()
    hcounts = hvals.shape[1]  # hmeta_.loc[hratio.columns, 'TISSUE'].value_counts()
    xlab = 'GTEx log2(Y/X expression)\nmean of n={} samples'.format(gcounts)
    ylab = 'HPA log2(Y/X expression)\nmean of n={} samples'.format(hcounts)
    
    if show_hpa_points:
        for p in _df.index:
            x = _df.at[p, 'gtex']
            ys = hvals.loc[p].dropna()
            y0 = np.min(ys)
            y1 = np.max(ys)
            ax.plot([x, x], [y0, y1], ls='-', lw=0.5, color='k', zorder=1)
        ylab += ' [min,max]'
    elif show_hpa_ci:
        for p in _df.index:
            x = _df.at[p, 'gtex']
            _yvals = hvals.loc[p].dropna()
            y0, y1 = np.percentile(_yvals, (5, 95))
            ax.plot([x, x], [y0, y1], ls='-', lw=0.5, color='k', zorder=1)
        ylab += ' [5th, 95th percentiles]'
    if show_gtex_ci:
        for p in _df.index:
            y = _df.at[p, 'hpa']
            _xvals = gvals.loc[p].dropna()
            x0, x1 = np.percentile(_xvals, (5, 95))
            ax.plot([x0, x1], [y, y], ls='-', lw=0.5, color='k', zorder=1)
        xlab += ' [5th, 95th percentiles]'
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    lim = pf.square_axes(ax)
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.plot(lim, lim, ls='--', color='0.5', zorder=1, lw=lw)

    r, p = ss.pearsonr(_df.gtex, _df.hpa)
    _ = ax.set_title('{0}: r2={1:.3f}, p={2:.2e}'.format(tissue, r**2, p))
    
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)
    for side in ['left', 'bottom']:
        ax.spines[side].set_linewidth(0.5)
    ax.tick_params(width=0.5)
    
    if smallfig:
        fig.set_figwidth(2.5)
        fig.set_figheight(2.5)
        ax = pf.format_axis_labels(ax, ytlabsize=7, titlesize=8)
    
    return fig, ax