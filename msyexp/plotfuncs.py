""" functions for plotting """
import numpy as np
import pandas as pd
import scipy.stats as ss
import matplotlib.pyplot as plt
import seaborn as sns

TISSUE_ABBR = {
    'Adipose - Subcutaneous' : 'Adipose (Subc)',
    'Adipose - Visceral'     : 'Adipose (Visc)',
    'Adrenal Gland'          : 'Adrenal Gland',
    'Artery - Aorta'         : 'Artery (Aort)',
    'Artery - Coronary'      : 'Artery (Coro)',
    'Artery - Tibial'        : 'Artery (Tib)',
    'Brain - Amygdala'       : 'Brain (Amyg)',
    'Brain - Cerebellum'     : 'Brain (Cblm)',
    'Brain - Cortex'         : 'Brain (Cort)',
    'Brain - Hippocampus'    : 'Brain (Hipp)',
    'Brain - Hypothalamus'   : 'Brain (Hypo)',
    'Brain - Striatum'       : 'Brain (Stri)',
    'Brain - Substantia nigra':'Brain (Subn)',
    'Breast'                 : 'Breast',
    'Colon - Sigmoid'        : 'Colon (Sigm)',
    'Colon - Transverse'     : 'Colon (Trns)',
    'Esophagus - Mucosa'     : 'Esoph (Muco)',
    'Esophagus - Muscularis' : 'Esoph (Musc)',
    'Heart - Atrial Appendage':'Heart (AtrA)',
    'Heart - Left Ventricle' : 'Heart (LVen)',
    'Kidney':'Kidney', 'Liver':'Liver', 'Lung':'Lung',
    'Nerve':'Nerve', 'Pancreas':'Pancreas', 
    'Pituitary':'Pituitary', 'Prostate':'Prostate',
    'Salivary Gland'         : 'Saliv Gland',
    'Skeletal Muscle'        : 'Sk Muscle',
    'Skin'                   : 'Skin',
    'Small Intestine'        : 'Sm Intestine',
    'Spinal Cord'            : 'Spinal Cord',
    'Spleen':'Spleen', 'Stomach':'Stomach',
    'Testis':'Testis', 'Thyroid':'Thyroid',
    'Cervix - Endocervix':'Cervix', 'Cervix':'Cervix',
    'Fallopian Tube'         :'Fallopian T',
    'Ovary':'Ovary', 'Uterus':'Uterus', 'Vagina':'Vagina'}

TISSUE_ORDER = ['Thyroid',
 'Kidney',
 'Prostate',
 'Bladder',
 'Salivary Gland',
 'Small Intestine',
 'Colon - Transverse',
 'Stomach',
 'Skin',
 'Esophagus - Mucosa',
 'Pancreas',
 'Adipose - Visceral',
 'Adipose - Subcutaneous',
 'Breast',
 'Lung',
 'Nerve',
 'Artery - Coronary',
 'Artery - Aorta',
 'Artery - Tibial',
 'Esophagus - Muscularis',
 'Colon - Sigmoid',
 'Adrenal Gland',
 'Heart - Left Ventricle',
 'Heart - Atrial Appendage',
 'Skeletal Muscle',
 'Spleen',
 'Liver',
 'Brain - Hippocampus',
 'Brain - Amygdala',
 'Brain - Striatum',
 'Brain - Hypothalamus',
 'Brain - Substantia nigra',
 'Brain - Cortex',
 'Spinal Cord',
 'Brain - Cerebellum',
 'Pituitary',
 'Testis']

TISSUE_ORDER_MF = ['Salivary Gland',
 'Breast',
 'Prostate',
 'Thyroid',
 'Lung',
 'Adipose - Visceral',
 'Adipose - Subcutaneous',
 'Nerve',
 'Esophagus - Muscularis',
 'Colon - Sigmoid',
 'Uterus',
 'Artery - Coronary',
 'Artery - Aorta',
 'Artery - Tibial',
 'Ovary',
 'Adrenal Gland',
 'Heart - Left Ventricle',
 'Heart - Atrial Appendage',
 'Skeletal Muscle',
 'Spleen',
 'Vagina',
 'Esophagus - Mucosa',
 'Skin',
 'Small Intestine',
 'Colon - Transverse',
 'Stomach',
 'Kidney',
 'Pancreas',
 'Liver',
 'Brain - Hippocampus',
 'Brain - Amygdala',
 'Brain - Cortex',
 'Brain - Hypothalamus',
 'Brain - Striatum',
 'Spinal Cord',
 'Brain - Substantia nigra',
 'Brain - Cerebellum',
 'Pituitary',
 'Testis']

# this list separates female and male-specific tissues
TISSUE_ORDER_MF2 = ['Uterus',
 'Ovary',
 'Vagina',
 'Salivary Gland',
 'Breast',
 'Thyroid',
 'Lung',
 'Adipose - Visceral',
 'Adipose - Subcutaneous',
 'Nerve',
 'Esophagus - Muscularis',
 'Colon - Sigmoid',
 'Artery - Coronary',
 'Artery - Aorta',
 'Artery - Tibial',
 'Adrenal Gland',
 'Heart - Left Ventricle',
 'Heart - Atrial Appendage',
 'Skeletal Muscle',
 'Spleen',
 'Esophagus - Mucosa',
 'Skin',
 'Small Intestine',
 'Colon - Transverse',
 'Stomach',
 'Kidney',
 'Pancreas',
 'Liver',
 'Brain - Hippocampus',
 'Brain - Amygdala',
 'Brain - Cortex',
 'Brain - Hypothalamus',
 'Brain - Striatum',
 'Spinal Cord',
 'Brain - Substantia nigra',
 'Brain - Cerebellum',
 'Pituitary',
 'Prostate',
 'Testis']

# sequential color palette
COLPAL_SEQ = sns.cubehelix_palette(50, light=0.95)

COLOR_BLUE = '#0570b0'
COLOR_DARKBLUE = '#023858'
COLOR_MIDNIGHT = '#10688d'
COLOR_COPPER = '#da7b28'
COLOR_BROWN = '#63300d'

def format_spines(ax, bottom=True):
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)
    for side in ['bottom', 'left']:
        ax.spines[side].set_color('k')
        ax.spines[side].set_linewidth(0.5)
    if not bottom:
        ax.spines['bottom'].set_visible(False)
        ax.tick_params(bottom=False)
    ax.tick_params(color='k', length=2, width=0.5)
    return ax

def format_axis_labels(ax, labsize=7, xtlabsize=6, ytlabsize=6, 
                       ticklabsize=None, titlesize=7, color='k'):
    """ set the size of tick and axis labels and color """
    if ticklabsize is not None:
        xtlabsize = ytlabsize = ticklabsize
    for lab in ax.get_yticklabels():
        lab.set_size(ytlabsize)
        lab.set_color(color)
    for lab in ax.get_xticklabels():
        lab.set_size(xtlabsize)
        lab.set_color(color)
    ax.set_ylabel(ax.get_ylabel(), size=labsize, color=color)
    ax.set_xlabel(ax.get_xlabel(), size=labsize, color=color)
    ax.set_title(ax.get_title(), size=titlesize, color=color)
    ax.tick_params(axis='both', pad=2)
    return ax

def simplify_tissues(ax, axis):
    """ Simplify tissue names in tick labels """
    tdict = TISSUE_ABBR
    if not axis in ('x', 'y'):
        raise ValueError("axis must be x or y")
    if axis == 'x':
        tlabs = [tdict[lab.get_text()] for lab in ax.get_xticklabels()]
        ax.set_xticklabels(tlabs)
    else:
        tlabs = [tdict[lab.get_text()] for lab in ax.get_yticklabels()]
        ax.set_yticklabels(tlabs)
    return ax

def rotate_ticklabels(ax, axis='x', rotation='vertical'):
    if axis == 'x':
        for lab in ax.get_xticklabels():
            lab.set_rotation(rotation)
    elif axis == 'y':
        for lab in ax.get_yticklabels():
            lab.set_rotation(rotation)
    else:
        raise ValueError("axis must be 'y' or 'x'")
    return ax

def square_axes(ax, equal_ticks=True, tick_axis='min'):
    """ Force x-axis and y-axis limits to be the same

    equal_ticks : bool
        require x- and y-axes to have the same ticks
    tick_axis : {'min'|'max'|'x'|'y'}
        how to select the ticks to enforce on both axes

    Returns new axis limits
    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    lim = (np.min((xlim[0], ylim[0])), np.max((xlim[1], ylim[1])))
    ax.set_xlim(lim)
    ax.set_ylim(lim)

    if equal_ticks:
        yticks = ax.get_yticks()
        xticks = ax.get_xticks()
        if tick_axis == 'min':
            ticks = yticks if len(yticks) < len(xticks) else xticks
        elif tick_axis == 'max':
            ticks = yticks if len(yticks) > len(xticks) else xticks
        elif tick_axis == 'x':
            ticks = xticks
        elif tick_axis == 'y':
            ticks = yticks
        else:
            raise ValueError("tick_axis must be 'min', 'max', 'x', or 'y'")
        ax.set_xticks(ticks)
        ax.set_yticks(ticks)

    return lim

def plot_male_female_stripplot(data, meta, gene, lo_thresh=0.1, log=True, 
                               palette={'male':COLOR_BLUE, 'female':'0.6'},
                               s=2.5, alpha=0.4):
    """ Plot stripplot of gene's expression levels in male/female samples
    """
    lab = '{}_TPM'.format(gene)
    df = pd.DataFrame({lab:data.loc[gene]})
    df[df < lo_thresh] = lo_thresh
    df = pd.merge(df, meta[['SAMPID', 'SEX', 'TISSUE']], left_index=True, 
                  right_on='SAMPID')

    fig, ax = plt.subplots(figsize=(6.5, 3))
    fig.subplots_adjust(bottom=0.28, top=0.92, left=0.08, right=0.85)

    ax = sns.stripplot(x='TISSUE', y=lab, hue='SEX', data=df,
                       dodge=True, alpha=alpha, s=s, 
                       palette=palette, order=TISSUE_ORDER_MF2, 
                       hue_order=['female', 'male'])

    if log:
        ax.set_yscale('log')

    ax = rotate_ticklabels(ax)
    ax = simplify_tissues(ax, 'x')
    ax = format_axis_labels(ax, ticklabsize=7)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
    ax.grid(axis='y', ls=':', lw=0.5, color='0.7')
    ax = format_spines(ax, bottom=False)
    ax.set_xlabel('')

    return fig, ax

def plot_XY_stripplot(data, meta, xgene, ygene, lo_thresh=0.1, log=True, 
                      xcolor=COLOR_COPPER, ycolor=COLOR_MIDNIGHT,
                      s=2.25, alpha=0.4, show_medians=True,
                      xcolor_mline=COLOR_BROWN, ycolor_mline=COLOR_DARKBLUE,
                      mline_lw=1):
    """ Plot stripplot of X vs. Y homolog expression in XY samples
    """
    df = data.loc[[xgene, ygene]].copy()
    df[df < lo_thresh] = lo_thresh
    df.index.name = 'gene'
    df = pd.melt(df.reset_index(), id_vars='gene', var_name='SAMPID', 
                 value_name='TPM')
    df = pd.merge(df, meta[['SAMPID', 'SEX', 'TISSUE']], on='SAMPID')

    fig, ax = plt.subplots(figsize=(6.5, 3))
    fig.subplots_adjust(bottom=0.28, top=0.92, left=0.08, right=0.85)

    tissues = df.TISSUE.unique()
    t_order = list(filter(lambda t: t in tissues, TISSUE_ORDER))
    ax = sns.stripplot(x='TISSUE', y='TPM', hue='gene', data=df,
                       dodge=True, alpha=alpha, s=s, 
                       palette={xgene:xcolor, ygene:ycolor}, 
                       order=t_order, 
                       hue_order=[xgene, ygene])
    if show_medians:
        mline_len = 0.4
        xmed = df.loc[df.gene==xgene].groupby('TISSUE').median()['TPM']
        ymed = df.loc[df.gene==ygene].groupby('TISSUE').median()['TPM']

        xlim = ax.get_xlim()
        ax.set_xlim(xlim)

        for i, t in enumerate(t_order):
            ax.plot([i-mline_len, i], [xmed[t], xmed[t]], 
                    lw=mline_lw, color=xcolor_mline, zorder=3)
            ax.plot([i+mline_len, i], [ymed[t], ymed[t]],
                    lw=mline_lw, color=ycolor_mline, zorder=3)

    if log:
        ax.set_yscale('log')

    ax = rotate_ticklabels(ax)
    ax = simplify_tissues(ax, 'x')
    ax = format_axis_labels(ax, ticklabsize=7)
    ax.legend(loc='upper left', bbox_to_anchor=(1, 1), frameon=False)
    ax.grid(axis='y', ls=':', lw=0.5, color='0.7')
    ax = format_spines(ax, bottom=False)
    ax.set_xlabel('')

    return fig, ax

def heatmap_2d(data1, data2, palette=None, vmin=None, vmax=None, 
    split_type='triangle', pad=0.05, pal_bins=50, na_color='0.7', 
    fill_na=False, ax=None, ax_cbar=None, plot_cbar=True, corner='upper left',
    **kwargs):
    """ Plot heatmap with two items per cell.
    
    Parameters
    ----------
    data1, data2 : pd.DataFrame
        dataframes of numerical values, with identical sets of columns and
        corresponding rows
    thresh : float
        cells exceeding this value will be plotted with a color
    palette : str
        name of matplotlib/seaborn color palette to use; if none given,
        will use seaborn default; if 'cubehelix', use cubehelix palette
    vmin, vmax : float
        minimum and maximum values for scaling colors 
    split_type : 'triangle'|'upper_lower'
        manner in which to divide each cell
    pad : 0 <= float < 0.5
        padding between boxes in heatmap
    pal_bins : int
        number of divisions of color palette
    na_color : colorspec
        color to assign to null values
    fill_na : bool
        if True, fill a NaN cell from data[1|2] with the cell from the other
        frame, if possible
    ax : Axes_Subplot
        existing Axes object to plot heatmap on
    ax_cbar : Axes_Suplot
        existing Axes object to plot color bar on
    plot_cbar : bool
        whether to plot color bar; this will be ignored if object is passed
        to ax argument but not ax_cbar argument
    corner : {'upper left'|'upper right'|'lower left'|'lower right'}
        when split_type is 'triangle', defines the corner of each cell in
        which the value from data1 will be plotted
    **kwargs : dict
        keyword arguments to pass to sns.cubehelix_palette for making a 
        custom palette
    
    Returns fig, (ax, ax_cbar)
    """
    if data1.shape != data2.shape:
        raise ValueError("dimensions of data1, data2 must be the same")
    if set(data1.columns) != set(data2.columns):
        raise ValueError("data1, data2 must have matching columns")
    if pad < 0 or pad >= 0.5:
        raise ValueError("pad must fall in range [0, 0.5)")
    if ax is None and ax_cbar is None:
        fig = plt.figure()
        if plot_cbar:
            ax = fig.add_axes((0.25, 0.25, 0.63, 0.65))
            ax_cbar = fig.add_axes((0.9, 0.6, 0.04, 0.3))
        else:
            ax = fig.add_subplot(111)
            ax_cbar = None
    elif ax is not None:
        fig = ax.get_figure()
    nrow, ncol = data1.shape
    ax.set_xlim((0, ncol))
    ax.set_ylim((0, nrow))
    wh = 1 - 2 * pad
    
    # reverse data1/2 index order for plotting convenience
    data1 = data1.loc[data1.index[::-1]]
    data2 = data2.loc[data2.index[::-1]]
    
    # set up color palette
    allvals = np.array(list(data1.values.flatten()) + list(data2.values.flatten()))
    allvals = allvals[~np.isnan(allvals)]
    if (allvals < 0).any():
        negvals = True
    else:
        negvals = False
    if palette is None:
        if negvals:
            pal = sns.color_palette('RdBu_r', pal_bins)
        else:
            pal = sns.cubehelix_palette(pal_bins)
    elif isinstance(palette, str):
        if palette == 'cubehelix':
            pal = sns.cubehelix_palette(pal_bins, **kwargs)
        else:
            pal = sns.color_palette(palette, pal_bins)
    else:
        pal = palette
    # duplicate last value of palette for mapping purposes
    pal = pal + [pal[-1]]
    # get min and max values for palette
    if (vmin is None) or (vmax is None):
        # set at least one of vmin, vmax automatically
        if vmin is not None:
            # need to set vmax only
            vmax = np.max(allvals)
        elif vmax is not None:
            # need to set vmin only
            vmin = np.min(allvals)
        else:
            # set both
            if negvals:
                vmax = np.max(np.abs(allvals))
                vmin = -1 * vmax
            else:
                vmax = np.max(allvals)
                vmin = np.min(allvals)
    vdiff = float(vmax - vmin)
    # truncate values
    data1[data1 < vmin] = vmin
    data2[data2 < vmin] = vmin
    data1[data1 > vmax] = vmax
    data2[data2 > vmax] = vmax
    # define function for mapping values to color
    def _map_color(x, vmin=vmin, vmax=vdiff, pal=pal, bins=pal_bins, na_color=na_color):
        return pal[int((x - vmin) / vdiff * bins)]
    
    if ax_cbar is not None:
        # plot colorbar
        cvals = np.linspace(vmin, vmax, pal_bins)
        cdiff = cvals[1] - cvals[0]
        for cv in cvals:
            p = plt.Rectangle([0, cv], 1, cdiff, edgecolor='none',
                              facecolor=_map_color(cv))
            ax_cbar.add_patch(p)
        ax_cbar.set_xlim((0, 1))
        ax_cbar.set_ylim((0, cvals[-1] + cdiff))
        ax_cbar.yaxis.tick_right()
        ax_cbar.tick_params(length=0)
        ax_cbar.set_yticks(np.linspace(vmin, vmax, 4))
        for side in ['top', 'bottom', 'left', 'right']:
            ax_cbar.spines[side].set_visible(False)
        ax_cbar.set_xticks([])
    
    for i in range(nrow):
        # cell vertex indices clockwise from bottom left
        y0 = y3 = i + pad
        y1 = y2 = i + pad + wh
        for j in range(ncol):
            x0 = x1 = j + pad
            x2 = x3 = j + pad + wh
            if np.isnan(data1.iat[i, j]):
                color1 = na_color
                isnan1 = True
            else:
                color1 = _map_color(data1.iat[i, j])
                isnan1 = False
            if np.isnan(data2.iat[i, j]):
                color2 = na_color
                isnan2 = True
            else:
                color2 = _map_color(data2.iat[i, j])
                isnan2 = False
            if fill_na and ((isnan1 + isnan2) == 1):
                # draw an unsplit square using the non-null value
                if isnan2:
                    color1 = _map_color(data1.iat[i, j])
                else:
                    color1 = _map_color(data2.iat[i, j])
                p1 = plt.Rectangle([x0, y0], wh, wh, edgecolor='none',
                                   facecolor=color1)
                ax.add_patch(p1)
            else:
                # create patches (in future could try other shapes)
                if split_type == 'triangle':
                    v0 = [x0, y0]
                    v1 = [x1, y1]
                    v2 = [x2, y2]
                    v3 = [x3, y3]
                    if corner == 'upper left':
                        coords1 = np.array([v0, v1, v2])
                        coords2 = np.array([v0, v2, v3])
                        # coords1 = np.array([[x0, y0], [x1, y1], [x2, y2]])
                        # coords2 = np.array([[x0, y0], [x2, y2], [x3, y3]])
                    elif corner == 'upper right':
                        coords1 = np.array([v1, v2, v3])
                        coords2 = np.array([v0, v1, v3])
                    elif corner == 'lower left':
                        coords1 = np.array([v0, v1, v3])
                        coords2 = np.array([v1, v2, v3])
                    elif corner == 'lower right':
                        coords1 = np.array([v0, v2, v3])
                        coords2 = np.array([v0, v1, v2])
                    else:
                        _opts = "'upper left', 'upper right', 'lower left', 'lower right'"
                        raise ValueError("corner argument must be one of %s" % _opts)
                    p1 = plt.Polygon(coords1, edgecolor='none', facecolor=color1)
                    p2 = plt.Polygon(coords2, edgecolor='none', facecolor=color2)
                elif split_type == 'upper_lower':
                    p1 = plt.Rectangle([x0, y0 + 0.5 * wh], wh, 0.5 * wh, 
                                       edgecolor='none', facecolor=color1)
                    p2 = plt.Rectangle([x0, y0], wh, 0.5 * wh,
                                       edgecolor='none', facecolor=color2)
                else:
                    raise ValueError("unrecognized value for split_type '%s'" % split_type)            
                ax.add_patch(p1)
                ax.add_patch(p2)

    ax.set_xticks(np.arange(ncol) + 0.5)
    ax.set_xticklabels(data1.columns, rotation='vertical')
    
    ax.set_yticks(np.arange(nrow) + 0.5)
    if corner == 'upper left':
        idx_str = "{0} / {1}"
    elif corner == "upper right":
        idx_str = "{1} \\ {0}"
    elif corner == "lower left":
        idx_str = "{0} \\ {1}"
    else:
        idx_str = "{1} / {0}"
    ytlabs = [idx_str.format(l1, l2) for l1, l2 in zip(data1.index, data2.index)]
    ax.set_yticklabels(ytlabs, style='italic')

    for side in ('top', 'right', 'bottom', 'left'):
        ax.spines[side].set_visible(False)
    
    return fig, (ax, ax_cbar)
