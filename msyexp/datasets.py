""" Functions for loading data """
import pandas as pd

from . import paths
from . import xytools
from . import general

def _load_dataframe(fp):
    if fp.endswith('h5'):
        data = pd.read_hdf(fp, 'df')
    else:
        data = pd.read_csv(fp, index_col=0, sep='\t')
    return data


def get_gtex_data(males_only=True, adjusted=False, rnaseqc=False,
                  sum_mc_families=True, combine_donor_reps=True,
                  transcripts=False):
    """ Return dataframe (genes, samples) of expression level estimates


    Parameters
    ----------
    males_only : bool
        if True, return only samples from male donors
    rnaseqc : bool
        if True, return expression level estimates reported by the 
        GTEx Consortium in their analysis with RNA-SeQC
    adjusted : bool
        if True, return covariate-adjusted data
    sum_mc_families : bool
        if True, sum the expression levels of multi-copy MSY genes
        by family
    combine_donor_reps : bool
        if True, when multiple samples are present from the same donor and
        tissue, avereage their expression levels
    combine_xy : bool
        return dataframe in which the expression levels of single-copy
        X- and Y-linked genes have been summed
    transcripts : bool
        return transcript-level TPMs
    """
    if adjusted:
        if rnaseqc:
            raise ValueError("no adjusted RNA-SeQC data available")
        if transcripts:
            raise ValueError("no adjusted transcript data available")

    # load dataset
    if rnaseqc:
        data = _load_dataframe(paths.EXP_RNASEQC)
    else:
        if transcripts:
            data = _load_dataframe(paths.EXP_TX_KSTO)
        else:
            if adjusted:
                if males_only:
                    data = _load_dataframe(paths.EXP_KSTO_ADJ_M)
                else:
                    data = _load_dataframe(paths.EXP_KSTO_ADJ_MF)
            else:
                data = _load_dataframe(paths.EXP_KSTO)

    meta = get_metadata()
    meta_ = meta.set_index('SAMPID')
    data = data.loc[:, meta_.loc[data.columns, 'TISSUE'].notnull()]
    excl_tissues = ['Bladder', 'Fallopian Tube', 'Cervix']
    data = data.loc[:, ~meta_.loc[data.columns, 'TISSUE'].isin(excl_tissues)]

    if males_only:
        data = data.loc[:, meta_.loc[data.columns, 'SEX']=='male']
    if combine_donor_reps:
        data = general.combine_donor_reps(data, meta)
    if (not transcripts) and sum_mc_families:
        data = xytools.sum_mc_families(data)

    return data


def get_metadata(study='gtex'):
    """ Return dataframe of metadata 

    Parameters
    ----------
    study : str {'gtex'|'hpa'}
        whether to return GTEx or HPA metadata
    """
    if study == 'gtex':
        meta = pd.read_csv(paths.METADATA_GTEX, low_memory=False, 
                           sep='\t')
    elif study == 'hpa':
        meta = pd.read_csv(paths.METADATA_HPA, sep='\t')
    else:
        raise ValueError("study must be one of 'gtex', 'hpa'")
    return meta


def get_anno_data(idx=None):
    """ Return dataframe with gene/transcript information """
    anno = pd.read_csv(paths.ANNOFILE, sep='\t')
    if idx is not None:
        anno = anno.drop_duplicates(idx).set_index(idx)
        if idx == 'gene_name':
            # add multicopy family names
            cols = [c for c in anno.columns]
            ygenes = filter(lambda g: not g in anno.index, xytools.XYMCFAM_NAMES)
            xgenes = filter(lambda g: not g in anno.index, ['VCX','HSFX'])
            df_y = pd.DataFrame({'seqname':'chrY', 'gene_type':'protein_coding'}, 
                                index=ygenes)
            df_x = pd.DataFrame({'seqname':'chrX', 'gene_type':'protein_coding'}, 
                                index=xgenes)
            df = pd.concat([df_y,df_x], axis=0, sort=False)
            anno = pd.concat([anno, df], axis=0, sort=False)
            anno = anno.loc[:, cols]
    return anno

def get_mappability():
    """ Return dataframe containing mappability by gene
    """
    mapp = pd.read_csv(paths.MAPP_DATA_V19, sep='\t')
    return mapp 

def get_simulation_results(sum_mc_families=True):
    """ Return dataframe with simulation results
    """
    data = {}
    sim_paths = [paths.SIM_Y0, paths.SIM_Y1_X2, 
                 paths.SIM_Y5_X10, paths.SIM_Y5_XUNCORR]
    labels = ['sim_y0', 'sim_y1_x2', 'sim_y5_x10', 'sim_y5_xuncorr']
    for path, sim in zip(sim_paths, labels):
        data[sim] = pd.read_csv(path, index_col=0, header=[0,1], sep='\t')
    data = pd.concat(data, axis=1)
    if sum_mc_families:
        data = xytools.sum_mc_families(data)
    return data

def get_hpa_data(males_only=True, sum_mc_families=True):
    """ Return dataframe (genes, samples) of HPA expression level estimates

    Parameters
    ----------
    males_only : bool
        if True, return only samples from male donors
    sum_mc_families : bool
        if True, sum the expression levels of multi-copy MSY genes
        by family
    """
    data = pd.read_csv(paths.EXP_HPA_KSTO, index_col=0, sep='\t')
    if males_only:
        meta = get_metadata('hpa')
        meta_ = meta.set_index("SAMPID")
        data = data.loc[:, meta_.loc[data.columns, 'SEX']=='male']
    if sum_mc_families:
        data = xytools.sum_mc_families(data)
    return data

def get_limma_sexbias():
    """ Return dataframe of results from limma sex-bias analysis """
    res = pd.read_csv(paths.LIMMA_SEXBIAS, sep='\t')
    return res

def get_mir1_expression():
    """ Return dataframe with miR-1 expression across human tissues
    from Ludwig et al. NAR (2016)
    """
    df = pd.read_csv(paths.MIR1_EXP, sep='\t')
    return df

def get_mass_spec_data():
    """ Return dataframes with peptide, protein, and metadata
    """
    dpep = pd.read_hdf(paths.MASS_SPEC_PEPTIDES, 'df')
    dprot = pd.read_hdf(paths.MASS_SPEC_PROTEINS, 'df')
    meta = pd.read_csv(paths.MASS_SPEC_METADATA, sep='\t')
    return dpep, dprot, meta

def get_bramerk_data():
    """ Return meta dataframe and Data dict with expression-level
    dataframes keyed by species

    Returns
    -------
    Data : dict, with keys {'chicken', 'chimp', 'mouse', 'rhesus'}
    meta : dataframe
    """
    meta = pd.read_csv(paths.BRAMERK_META, sep='\t')
    Data = {}
    spc_list = ('chicken', 'chimp', 'mouse', 'rhesus')
    path_list = (paths.BRAMERK_CHICKEN, paths.BRAMERK_CHIMP, 
             paths.BRAMERK_MOUSE, paths.BRAMERK_RHESUS)
    for spc, path in zip(spc_list, path_list):
        data_s = pd.read_csv(path, sep='\t', index_col=0)
        
        if spc == 'mouse':
            data_s = data_s.rename(index={'ENSMUST00000087143.6':'EIF1AX'})
        elif spc == 'chicken':
            data_s = data_s.rename(index={'ENSGALT00000026468.6':'EIF1AX'})
            
        Data[spc] = data_s
    return Data, meta