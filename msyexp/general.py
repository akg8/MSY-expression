""" General calculations/transformations """
import sys

import numpy as np 
import pandas as pd
import scipy.stats as ss  

def _print(msg, newline=True):
    sys.stdout.write(msg)
    if newline:
        sys.stdout.write('\n')
    sys.stdout.flush()

def split_male_female_samples(data, meta):
    """ Return dataframes with male and female data respectively
    """
    meta_ = meta.set_index('SAMPID')
    mdata = data.loc[:, meta_.loc[data.columns, 'SEX']=='male']
    fdata = data.loc[:, meta_.loc[data.columns, 'SEX']=='female']
    return mdata, fdata

def stat_by_tissue(data, meta, tissuecol='TISSUE', func='median', runcol='SAMPID'):
    """ Return summary statistics for genes in tissues. 
    
    Parameters
    ----------
    data : pandas.DataFrame or pandas.Series
      genes x samples dataframe with gene ID in index, or a series object
      of expression values indexed by sample ID
      
    meta : pandas.DataFrame
      preprocessed metadata dataframe; sample IDs should be in 'RUN' column,
      and donor IDs in 'SUBJID' column
    
    tissuecol : str
      specifies how biological replicates are defined; should match the
      name of a column in the metadata table
      
    func : str, ['mean'|'median'|'geomean']
      a summary statistic to calculate for each gene in each tissue
      
    Returns
    -------
    res : pandas.DataFrame or pandas.Series
      If data is a DataFrame, will return a genes x tissue dataframe
      containing the summary statistic for each gene in each tissue;
      if data is a Series, will return a series indexed by tissue
    """
    series = False
    if isinstance(data, pd.Series):
        series = True
        data = pd.DataFrame(data).T
    meta_ = meta.set_index(runcol)
    if func == 'median':
        res = data.groupby(by=lambda r: meta_.at[r, tissuecol], axis=1).median()
    elif func == 'mean':
        res = data.groupby(by=lambda r: meta_.at[r, tissuecol], axis=1).mean()
    elif func == 'geomean':
        # geometric mean
        data = np.log2(data + 0.5)
        res = data.groupby(by=lambda r: meta_.at[r, tissuecol], axis=1).mean()
        res = (2**res) - 0.5
    else:
        raise ValueError("unrecognized value for func argument '%s'" % func)
    if series:
        res = res.iloc[0]
    return res

def combine_donor_reps(data, meta, tissuecol='TISSUE', runcol='SAMPID',
                       tpm_renorm=True):
    """ Return dataframe with mean of donor replicates.
    
    For each set of donor replicates (e.g. samples from the same donor and
    tissue, the sample name given to the combined values in the returned 
    dataframe is equal to the first encountered sample from that set of 
    replicates.
    
    Parameters
    ----------
    data : pandas.DataFrame
      genes x samples dataframe with gene ID in index
      
    meta : pandas.DataFrame
      preprocessed metadata dataframe; sample IDs should be in 'RUN' column,
      and donor IDs in 'SUBJID' column
    
    tissuecol : str
      specifies how biological replicates are defined; should match the
      name of a column in the metadata table

    runcol : str
      column of meta containing sample IDs

    tpm_renorm : bool
      renormalize expression levels to TPM units
    """
    meta_ = meta.set_index(runcol)
    # sort data columns for consistency
    cols = [c for c in data.columns]
    cols.sort()
    data = data.loc[:, cols]
    
    # create donor-tissue IDs
    run2id = {}
    id2run = {}
    for run in data.columns:
        idx = "{0}_{1}".format(meta_.at[run, 'SUBJID'],
                               meta_.at[run, tissuecol])
        run2id[run] = idx
        if not idx in id2run:
            id2run[idx] = run
    
    # combine values
    data = data.groupby(by=lambda r: run2id[r], axis=1).mean()
    if tpm_renorm:
        data = data / data.sum(axis=0) * 1e6
    
    # reset column names
    data.columns = map(lambda c: id2run[c], data.columns)
    return data


def calc_hkeeping_norm_factor(data, exp_range=[10, 100], n_genes=50,
                              func='mean', ref=None, groups=None):
    """ Calculate normalization factor based on housekeeping genes 
    
    select <n_genes> with median exp levels in exp_range with most
    conserved ranks
    """
    data = data.copy()
    if ref is not None:
        if isinstance(ref, pd.Series):
            data['_ref'] = ref
            refcol = '_ref'
        else:
            refcol = ref
    else:
        refcol = None
    
    ## weighted by tissue
    if groups is not None:
        if func == 'mean':
            mbt = np.log2(data+0.5).groupby(by=lambda s: groups[s], axis=1).mean()
            meds = 2**(mbt.mean(axis=1)) - 0.5
        else:
            mbt = data.groupby(by=lambda s: groups[s], axis=1).median()
            meds = mbt.median(axis=1)
        
        # find genes in expression range
        genes = meds.loc[meds.between(exp_range[0], exp_range[1])].index
        if len(genes) < n_genes:
            raise ValueError("not enough genes in expression level range")
        
        # calculate ranks and noise based on tissue-level estimates
        ranks = mbt.rank(axis=0)
        rnoise = ranks.std(axis=1) / ranks.mean(axis=1)
            
    ## unweighted
    else:        
        if func == 'mean':
            meds = 2**((np.log2(data + 0.1)).mean(axis=1))
        else:
            meds = data.median(axis=1)
        
        # find genes in expression range
        genes = meds.loc[meds.between(exp_range[0], exp_range[1])].index
        if len(genes) < n_genes:
            raise ValueError("not enough genes in expression level range")
            
        # calculate ranks and noise based on full data
        ranks = data.rank(axis=0)
        rnoise = ranks.std(axis=1) / ranks.mean(axis=1)
    
    # select housekeeping genes
    genes = rnoise.loc[genes].sort_values()[:n_genes].index
    
    # calculate scale factors
    logfc = np.log(data.loc[genes].divide(meds.loc[genes], axis=0))
    scale = np.exp(logfc.median(axis=0))
    if refcol is not None:
        scale = scale / scale[refcol]
        if refcol == '_ref':
            _ = scale.pop('_ref')
    return scale, list(genes)


def housekeeping_normalize(data, exp_range=[10, 100], n_genes=50, ref=None,
                           weighted=False, meta=None):
    """ Normalize expression levels based on housekeeping genes

    Housekeeping genes are the <n_genes> genes with mean expression
    levels falling in the specified expression-level range that 
    have the most conserved ranks across the samples (ranks with lowest
    coefficient of variation [stdev / mean]).

    Parameters
    ----------
    data : pd.DataFrame
        genes x samples/tissues dataframe
    exp_range : list-like
        consider genes with mean expression levels that are 
        greater than exp_range[0] and less than exp_range[1]
    n_genes : int
        number of housekeeping genes used to calculate normalization factor
    ref : pd.Series or str
        the name of a column in data or a Series object (with the same index
        as data) to normalize to--after scale factors are calculated, they
        will be adjusted so that the scale for this sample is 1 
    weighted : bool
        if True, housekeeping genes are selected as those genes with the
        most conserved expression-level ranks across *tissues*, by giving
        more weight to samples collected from tissues where fewer samples
        were collected; otherwise, all samples will be weighted equally,
        irrespective of tissue
    meta : pd.DataFrame
        metadata dataframe; needed for weighted calculation
    """
    if ref is not None:
        if isinstance(ref, str):
            if not ref in data.columns:
                raise ValueError("reference sample not found in columns")
        elif isinstance(ref, pd.Series):
            if set(ref.index) != set(data.index):
                raise ValueError("ref and data index do not match")
    if weighted:
        if meta is None:
            raise ValueError("must supply meta dataframe to determine weights")
        else:
            meta_ = meta.set_index('SAMPID')
            groups = meta_.loc[data.columns, 'TISSUE']
    else:
        groups = None
    
    hfac, genes = calc_hkeeping_norm_factor(data, exp_range=exp_range, 
                                            n_genes=n_genes, ref=ref, 
                                            groups=groups)
    data = data / hfac
    return data


def mannwhitney_permute(x1, x2, n_boot=10000):
    """ Two-sided Mann-Whitney U test by permutation
    """
    # calc U statistic
    n1 = len(x1)
    n2 = len(x2)
    mid = 0.5 * n1 * (n1 + 1)
    x = list(x1)
    x.extend(x2)
    x = np.array(x)
    if len(ss.find_repeats(x).values) > 0:
        raise ValueError("cannot compute exact p-value with ties")
    ranks = ss.rankdata(x)
    r1 = np.sum(ranks[:n1])
    u1 = r1 - mid
    u = np.min((u1, n1*n2 - u1))

    # permute
    u_dist = np.empty(n_boot)
    for i in range(n_boot):
        np.random.shuffle(x)
        ranks = ss.rankdata(x)
        r1 = np.sum(ranks[:n1])
        u1 = r1 - mid
        u_dist[i] = np.min((u1, n1*n2 - u1))
    pval = np.sum(u_dist <= u) / float(n_boot)
    return u, pval 
