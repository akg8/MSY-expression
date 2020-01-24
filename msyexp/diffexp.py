""" functions for co-expression and differential expression """
import numpy as np
import pandas as pd
import scipy.stats as ss
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns 
import statsmodels.api as sm 

from . import datasets as ds 
from . import general
from . import xytools
from . import plotfuncs as pf

### AUXILIARY FUNCTIONS ###

def residualize(data, covs):
	""" Residualize data using covariates 

	data : DataFrame (n_genes, n_indivs)
	covs : DataFrame (n_indivs, n_covars)
	"""
	# calculate hat matrix
	c = np.matrix(covs)
	hat = c*np.linalg.inv(np.transpose(c)*c)*np.transpose(c)

	X = np.matrix(data.T)
	X_ = X - hat * np.matrix(X)
	X_ = pd.DataFrame(X_, index=data.columns, columns=data.index).T
	return X_

def get_expression_PCs(data, covs=None, n_PCs=None, median_thresh=5.,
					   hkeep_norm=True):
	""" Calculate expression PCs after regressing out covariates 
	
	Parameters
	----------
	data : DataFrame
		genes x samples;
	covs : DataFrame
		explicitly modeled covariates (samples x variables); if not passed,
		will just subtract off the mean
	n_PCs : int
		number of PCs to be returned; if 'None', return all
	median_thresh : float
		only calculate the PCs using genes with a median TPM above this level;
		to skip this, set to a value less than 0
	hkeep_norm : bool
		apply housekeeping normalization
	"""
	D = data.copy()
	# set up covariates
	if covs is not None:
		if D.shape[1] != covs.shape[0]:
			raise ValueError("dimension of `data` and `covs` not compatible")
		if set(D.columns) != set(covs.index):
			raise ValueError("samples in `data` and `covs` are not identical")
		
		X = covs.copy()
		X = X.reindex(D.columns, axis=0)
	
		# add constant to X
		if not (covs==1.).all(axis=0).any():
			X['const'] = 1.
	else:
		X = pd.DataFrame({'const':1.}, index=D.columns)
	
	# prepare data
	if median_thresh > 0:
		okgenes = D.median(axis=1) > median_thresh
	else:
		okgenes = pd.Series(True, index=D.index)
	if hkeep_norm:
		D = general.housekeeping_normalize(D)
	D = D.loc[okgenes]
	D = np.log2(D + 0.5)
	R = residualize(D, X)
	
	# PCA
	pca = PCA().fit(R)
	if n_PCs is None:
		n_PCs = pca.n_components_
	cols = ['PC{}'.format(i+1) for i in range(n_PCs)]
	comps = pd.DataFrame(pca.components_[:n_PCs].T, 
						 index=D.columns, columns=cols)
	return comps

def regress(Y, covs, var_name):
    X = covs.copy()
    X = pd.concat([pd.DataFrame({'const':1}, index=X.index), X], axis=1)
    Y = Y.reindex(X.index, axis=1)
    B = pd.DataFrame(np.linalg.pinv(X).dot(Y.T),
                     index=X.columns, columns=Y.index)
    resid = Y - X.dot(B).T
    sigma_h2 = (resid**2).sum(axis=1) / float(len(Y.columns))
    Q = pd.DataFrame(np.linalg.inv(X.T.dot(X)), index=X.columns,
                     columns=X.columns)
    b_se = np.sqrt(Q.at[var_name, var_name] * sigma_h2)
    # t-statistics
    ts = B.loc[var_name] / b_se.loc[B.columns]
    # p-values
    pvals = 2*pd.Series(ss.t.cdf(-1*np.abs(ts), df=X.shape[0]-X.shape[1]-1), ts.index)

    res = pd.DataFrame({'log2fc':B.loc[var_name], 't':ts, 'pval':pvals})
    res = res.loc[:, ['log2fc', 't', 'pval']]
    _rej, pval_adj, _, _ = sm.stats.multipletests(res.pval, method='fdr_bh')
    res['pval_adj'] = pval_adj
    res = res.sort_values('pval_adj')
    return res


### ANALYSIS FUNCTIONS ###

def analyze_coexp_onetissue(dt, xgenes, ygenes, method='spearman',
							median_thresh=5, n_PCs=1):
	""" Analyze co-expression of X-Y gene pairs in one tissue 

	Parameters
	----------
	dt : DataFrame
		genes x samples dataframe for one tissue
	xgenes, ygenes : list
		list of X homologs and corresponding Y homologs
	method : {'spearman'|'pearson'}
		type of correlation to calculate
	median_thresh : float
		only analyze genes that have a median expression level
		above this value
	n_PCs : int
		number of expression PCs to control for

	Returns
	-------
	res_t : DataFrame
		long-form dataframe with coexpression results
	dt_r : DataFrame
		genes x samples dataframe of residualized expression levels
		used to calcualte correlations
	"""
	if len(xgenes) != len(ygenes):
		raise ValueError("xgenes and ygenes must be the same length")
	# filter and normalize
	okgenes = dt.median(axis=1) > median_thresh
	dt = general.housekeeping_normalize(dt)
	dt = dt.loc[okgenes]

	# get expression PC(s)
	comps = get_expression_PCs(dt, n_PCs=n_PCs, median_thresh=-1, 
							   hkeep_norm=False)
	## ensure centered
	comps = comps - comps.mean(axis=0)
	## comps['const'] = 1.

	# residualize data (remove effect of PC)
	dt = np.log2(dt + 0.5)
	## subtract off mean
	means = dt.mean(axis=1)
	dt = dt.subtract(means, axis=0)
	##
	dt_r = residualize(dt, comps)
	## add mean back
	dt_r = dt_r.add(means, axis=0)

	# calculate gene-gene correlations
	if method == 'spearman':
		cor_t = pd.DataFrame(np.corrcoef(dt_r.rank(axis=1)), 
							 index=dt_r.index, 
							 columns=dt_r.index)
		mlab = 'spearman'
	else:
		cor_t = pd.DataFrame(np.corrcoef(dt_r), index=dt_r.index, 
							 columns=dt_r.index)
		mlab = 'pearson'

	gx_t = {}
	gy_t = {}
	r_xy_avg_t = {}
	rx_y_t = {}
	ry_x_t = {}
	pear_r_t = {}
	pear_p_t = {}
	for gx, gy in zip(xgenes, ygenes):
		if (gx in cor_t.index) and (gy in cor_t.index):
			lab = '{0}/{1}'.format(gx, gy)
			gx_t[lab] = gx
			gy_t[lab] = gy

			# pearson r
			if method == 'spearman':
				pear_r_t[lab], pear_p_t[lab] = ss.spearmanr(dt_r.loc[gx], 
															dt_r.loc[gy])
			else:
				pear_r_t[lab], pear_p_t[lab] = ss.pearsonr(dt_r.loc[gx], 
														   dt_r.loc[gy])

			# reciprocal ranks
			rx = cor_t.loc[gx].rank(ascending=False) - 1
			ry = cor_t.loc[gy].rank(ascending=False) - 1
			rx_y_t[lab] = rx[gy]
			ry_x_t[lab] = ry[gx]

			# average rank (quantile)
			r_xy_avg_t[lab] = 0.5 * (rx[gy] + ry[gx]) / (len(rx) - 1)

	res_t = pd.DataFrame({'X_homolog':gx_t, 'Y_homolog':gy_t, 
						  '{}_r'.format(mlab):pear_r_t, 
						  '{}_p'.format(mlab):pear_p_t,
						  'rank_X(Y)':rx_y_t, 'rank_Y(X)':ry_x_t,
						  'avg_rank':r_xy_avg_t})
	res_t = res_t.reindex(['X_homolog', 'Y_homolog',
						   '{}_r'.format(mlab), '{}_p'.format(mlab),
						   'rank_X(Y)', 'rank_Y(X)', 'avg_rank'], axis=1)
	
	res_t.index.name = 'pair'
	res_t = res_t.reset_index()
	return res_t, dt_r

def analyze_coexpression(data, meta, xgenes=xytools.XYPAIRS_X9,
						 ygenes=xytools.XYPAIRS_Y9, method='spearman',
						 median_thresh=3, n_PCs=1, min_tissue_samps=30,
						 verbose=True):
	""" Analyze co-expression of X-Y gene pairs in all tissues

	Parameters
	----------
	data : DataFrame
		genes x samples dataframe for all tissues
	meta : DataFrame
		samples x variables metadata dataframe
	xgenes, ygenes : list
		list of X homologs and corresponding Y homologs
	method : {'spearman'|'pearson'}
		type of correlation to calculate
	median_thresh : float
		only analyze genes that have a median expression level
		above this value
	n_PCs : int
		number of expression PCs to control for
	min_tissue_samps : int
		only analyze tissues where there at least this many samples

	Returns
	-------
	res : DataFrame
		long-form dataframe with coexpression results
	data_r : DataFrame
		genes x samples dataframe of residualized expression levels
		used to calcualte correlations. note that this dataframe 
		will contain null values in cases where a gene was not 
		sufficienly highly expressed the samples from that tissue
		to be analyzed
	"""
	meta_ = meta.set_index('SAMPID')
	tcounts = meta_.loc[data.columns, 'TISSUE'].value_counts()
	tissues = tcounts.index[tcounts >= min_tissue_samps]

	res = []
	data_r = []
	for i, t in enumerate(tissues):
		if verbose:
			general._print("{0} [{1}/{2}]...".format(t, i+1, len(tissues)))
		dt = data.loc[:, meta_.loc[data.columns, 'TISSUE']==t].copy()
		res_t, dt_r = analyze_coexp_onetissue(dt, xgenes, ygenes, method=method,
											  median_thresh=median_thresh, 
											  n_PCs=n_PCs)
		cols = [c for c in res_t.columns]
		cols_ = cols[:3] + ['tissue'] + cols[3:]
		res_t['tissue'] = t
		res_t = res_t.reindex(cols_, axis=1)
		res.append(res_t)
		data_r.append(dt_r)
	res = pd.concat(res, axis=0, ignore_index=True)
	data_r = pd.concat(data_r, axis=1, sort=True)
	return res, data_r


def find_DEtissues_onegene(data, meta, gene, min_samps=30, 
						   min_l2fc=np.log2(1.5), p_thresh=1e-3, 
						   frac_tissues=0.75):
	""" Find tissues where a gene is <min_fc> up or downregulated
	relative to <frac_tissues> of other tissues at a p-value of 
	<p_thresh> 
	"""
	meta_ = meta.set_index('SAMPID')

	dg = np.log2(data.loc[gene] + 0.5)

	# subset to tissues with enough samples
	tcounts = meta_.loc[dg.index, 'TISSUE'].value_counts()
	tissues = list(tcounts.index[tcounts >= min_samps])
	tissues.sort()
	dg = dg.loc[meta_.loc[dg.index, 'TISSUE'].isin(tissues)]

	# get summary statistic
	dg_med = dg.groupby(by=lambda s: meta_.at[s, 'TISSUE']).median()

	# perform pairwise t-tests
	pvals = pd.DataFrame(np.zeros((len(tissues), len(tissues))),
						 index=tissues, columns=tissues)
	fcs = pd.DataFrame(np.zeros((len(tissues), len(tissues))),
					   index=tissues, columns=tissues)

	for i in range(len(tissues)):
		for j in range(i+1, len(tissues)):
			t1 = tissues[i]
			t2 = tissues[j]
			t_, p_ = ss.ttest_ind(dg.loc[meta_.loc[dg.index, 'TISSUE']==t1],
								  dg.loc[meta_.loc[dg.index, 'TISSUE']==t2],
								  equal_var=False)
			pvals.at[t1, t2] = pvals.at[t2, t1] = p_
			fcs.at[t1, t2] = dg_med.at[t1] - dg_med.at[t2]
			fcs.at[t2, t1] = -1 * fcs.at[t1, t2]
			
	for t in tissues:
		pvals.at[t, t] = np.nan
		fcs.at[t, t] = np.nan

	# collect results    
	res = {}
	_tup = pvals.lt(p_thresh) & fcs.gt(min_l2fc)
	_tdown = pvals.lt(p_thresh) & fcs.lt(-1*min_l2fc)
	res['n_up'] = _tup.sum(axis=1)
	res['n_down'] = _tdown.sum(axis=1)
	res['med_log2fc'] = fcs.median(axis=1)
	res = pd.DataFrame(res)
	sig_up = res.loc[:, 'n_up'] > ((len(tissues)-1)*frac_tissues)
	sig_dn = res.loc[:, 'n_down'] > ((len(tissues)-1)*frac_tissues)
	res['direction'] = 'n.s.'
	res.loc[sig_up, 'direction'] = 'up'
	res.loc[sig_dn, 'direction'] = 'down'

	return res

def find_DE_tissues_multi(data, meta, genes, min_samps=30, min_l2fc=np.log2(1.5),
						  p_thresh=1e-3, frac_tissues=0.75, add_chrom=True,
						  hkeep_norm=True, verbose=True):
	""" find tissues where a gene is <min_fc> up or downregulated
	relative to <frac_tissues> of other tissues at a p-value of 
	<p_thresh> 
	
	for multiple genes
	"""
	if hkeep_norm:
		data = general.housekeeping_normalize(data, weighted=True, meta=meta)

	res = []
	for i, g in enumerate(genes):
		if verbose:
			general._print("{0} [{1}/{2}]...".format(g, i+1, len(genes)))
		res_g = find_DEtissues_onegene(data, meta, g, min_samps=min_samps,
									   min_l2fc=min_l2fc, p_thresh=p_thresh,
									   frac_tissues=frac_tissues)
		res_g.index.name = 'tissue'
		res_g['gene'] = g
		res_g = res_g.reset_index()
		res.append(res_g)
	
	res = pd.concat(res, axis=0, ignore_index=True)
	cols = ['gene', 'tissue', 'med_log2fc', 'direction', 'n_up', 'n_down']
	if add_chrom:
		annog = ds.get_anno_data('gene_name')
		annog = annog.rename(columns={'seqname':'chrom'})
		res = pd.merge(res, annog[['chrom']], left_on='gene', 
					   right_index=True)
		cols = cols[:1] + ['chrom'] + cols[1:]
	res = res.reindex(cols, axis=1)

	# sort results
	res['fc_abs'] = np.abs(res.med_log2fc)
	res.loc[:, 'direction'] = pd.Categorical(res.loc[:, 'direction'], ['up', 'down', 'n.s.'])
	res = res.sort_values(['direction', 'fc_abs'], ascending=(True, False))
	_ = res.pop('fc_abs')
	
	return res

def plot_DE_results(deres, sort_by_pair=True):
    if sort_by_pair:
        genes = []
        for gx, gy in zip(xytools.XYPAIRS_X9, xytools.XYPAIRS_Y9):
            if gx in deres.gene.values and gy in deres.gene.values:
                genes.extend((gy, gx))
    else:
        genes = (deres.gene.unique())
        genes.sort()
    tissues = list(filter(lambda t: t in deres.tissue.values, pf.TISSUE_ORDER))
    
    # pivot log2fc values
    df = deres.loc[:, ['gene', 'tissue', 'med_log2fc']]
    df = df.pivot(index='gene', columns='tissue', values='med_log2fc')
    df = df.loc[genes, tissues]
    
    # pivot significance values
    df2 = deres.loc[:, ['gene', 'tissue', 'direction']]
    df2 = df2.pivot(index='gene', columns='tissue', values='direction')
    df2 = df2.loc[genes, tissues]
    
    fig, ax = plt.subplots(figsize=(7.5, 5))
    fig.subplots_adjust(bottom=0.25, right=1)
    ax = sns.heatmap(df, ax=ax, cmap='RdBu_r', center=0, 
                     cbar_kws={'label':'log2(FC) vs. median tissue'})
    ax.set_ylim(df.shape[0], 0)
    
    # plot dividing lines
    if sort_by_pair:
        for y_ in np.arange(2, len(genes), 2):
            ax.plot(ax.get_xlim(), (y_, y_), color='w', lw=1)
    
    for lab in ax.get_yticklabels():
        lab.set_style('italic')
    ax.tick_params(axis='both', left=False, bottom=False, pad=1)
    
    # mark significance
    xs = []
    ys = []
    for i in range(df2.shape[0]):
        xs_ = np.where(df2.iloc[i] != 'n.s.')[0]
        for x_ in xs_:
            ys.append(i)
            xs.append(x_)
    xs = np.array(xs) + 0.5
    ys = np.array(ys) + 0.5
    ax.scatter(xs, ys, edgecolors='none', facecolors='w', marker='*')
    ax.set_xlabel('')
    ax.set_ylabel('')

    ax = pf.simplify_tissues(ax, axis='x')
    
    return fig, ax

def plot_xx_vs_xy_stackedbar(data, meta, xgene, ygene, norm_xx=True, ax=None,
                             xcolor_m='0.7', xcolor_f='0.7', ycolor=pf.COLOR_BLUE,
                             bar_width=0.35, pad=0.1, err_quantiles = [0.25, 0.75],
                             show_xm_err=False, pvals=None, p_offset=0.15, p_size=9,
                             p_color='0.6'):
    """
    pvals : pd.Series
        Series of p-values equal in length to the nubmer of tissues
    """
    meta_ = meta.set_index('SAMPID')
    tissues = list(filter(lambda t: not t in ('Prostate', 'Testis', 'Kidney', 'Bladder'), pf.TISSUE_ORDER))
    if pvals is not None and not set(tissues).issubset(pvals.index):
        raise ValueError("some tissues not found in pvals index")

    df_x = data.loc[xgene]
    df_xy = data.loc[xgene] + data.loc[ygene]

    df_x = df_x.loc[meta_.loc[df_x.index, 'TISSUE'].isin(tissues)]
    df_xy = df_xy.reindex(df_x.index)

    xlab = '{}_TPM'.format(xgene)
    xylab = '{0}+{1}_TPM'.format(xgene, ygene)

    df = pd.DataFrame({xlab:df_x, xylab:df_xy})
    df = pd.merge(df, meta[['SAMPID', 'TISSUE', 'SEX']], left_index=True, right_on='SAMPID')

    if norm_xx:
        # calculate mean XX expression
        xx_mean = df.loc[df.SEX.eq('female'), [xlab, 'TISSUE']].groupby('TISSUE').mean()
        xx_mean.columns = ['XX_mean']
        df = pd.merge(df, xx_mean, left_on='TISSUE', right_index=True)
        df[xlab] = df[xlab] / df['XX_mean']
        df[xylab] = df[xylab] / df['XX_mean']

    # collect bar heights and errors
    # X in females
    xf_group = df.loc[df.SEX.eq('female'), [xlab, 'TISSUE']].groupby('TISSUE')
    x_f_m = xf_group.mean()[xlab].reindex(tissues)
    x_f_err0 = -1*(xf_group.quantile(err_quantiles[0])[xlab].reindex(tissues) - x_f_m)
    x_f_err1 = xf_group.quantile(err_quantiles[1])[xlab].reindex(tissues) - x_f_m

    # X in males
    xym_group = df.loc[df.SEX.eq('male'), [xlab, xylab, 'TISSUE']].groupby('TISSUE')
    xym_q0 = xym_group.quantile(err_quantiles[0])
    xym_q1 = xym_group.quantile(err_quantiles[1])
    x_m_m = xym_group.mean()[xlab].reindex(tissues)
    x_m_err0 = -1*(xym_q0[xlab].reindex(tissues) - x_m_m)
    x_m_err1 = xym_q1[xlab].reindex(tissues) - x_m_m

    # X+Y in males
    xy_m_m = xym_group.mean()[xylab].reindex(tissues)
    xy_m_err0 = -1*(xym_q0[xylab].reindex(tissues) - xy_m_m)
    xy_m_err1 = xym_q1[xylab].reindex(tissues) - xy_m_m
    
    if ax is None:
        fig, ax = plt.subplots(figsize=(9, 4))
    else:
        fig = ax.get_figure()

    xs = np.arange(len(tissues))
    error_kws = {'lw':0.5, 'color':0.1}
    # X in females
    ax.bar(xs - 0.5*(pad+bar_width), x_f_m, width=bar_width, 
           yerr=[x_f_err0.values, x_f_err1.values],
           error_kw=error_kws, edgecolor='none', facecolor=xcolor_f)

    # X+Y in males
    ax.bar(xs + 0.5*(pad+bar_width), xy_m_m, width=bar_width, 
           yerr=[xy_m_err0.values, xy_m_err1.values],
           error_kw=error_kws, edgecolor='none', facecolor=ycolor)

    # X in males
    if show_xm_err:
        ax.bar(xs + 0.5*(pad+bar_width), x_m_m, width=bar_width,
               yerr=[x_m_err0.values, x_m_err1.values],
               error_kw={'lw':0.5}, ecolor=xcolor_m,
               edgecolor='none', facecolor=xcolor_m)
    else:
        ax.bar(xs + 0.5*(pad+bar_width), x_m_m, width=bar_width,
               edgecolor='none', facecolor=xcolor_m)

    if pvals is not None:
        for i, t in enumerate(tissues):
            if pvals[t] < 0.05:
                if xy_m_m[t] > x_f_m[t]:
                    y_p = xy_m_m[t]
                else:
                    y_p = x_f_m[t]
                ax.text(i, y_p * (1 + p_offset), '*', size=p_size, 
                        color=p_color, horizontalalignment='center',
                        verticalalignment='center')

    ax.set_xticks(xs)
    xlim = [-0.75, len(tissues)-0.25]
    ax.set_xlim(xlim)
    ax.set_xticklabels(tissues, rotation='vertical')

    ax = pf.simplify_tissues(ax, 'x')
    ax = pf.format_spines(ax)

    if norm_xx:
        ax.plot(xlim, [0.5, 0.5], lw=0.5, ls=':', color='w')
        ax.plot(xlim, [1, 1], lw=0.5, ls=':', color='0.7', zorder=1)
        ax.plot(xlim, [2, 2], lw=0.5, ls=':', color='0.7', zorder=1)
        
    return fig, ax

