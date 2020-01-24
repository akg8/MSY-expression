""" Variables and functions for X and Y chromosome genes """

# MSY protein coding genes/gene families from Skaletsky et al. 2003
MSYFAMS = ['AMELY', 'BPY2', 'CDY', 'DAZ', 'DDX3Y', 'EIF1AY', 'HSFY',
           'KDM5D', 'NLGN4Y', 'PCDH11Y', 'PRKY',
           'PRY', 'RBMY', 
           'RPS4Y1', 'RPS4Y2', 'SRY', 'TBL1Y', 'TGIF2LY', 'TMSB4Y',
           'TSPY', 'TXLNGY', 'USP9Y', 'UTY', 'VCY', 'XKRY', 'ZFY']

# MSY protein-coding genes
MSYGENES = ['AMELY', 'BPY2', 'BPY2B', 'BPY2C', 'CDY1', 'CDY1B', 'CDY2A', 'CDY2B',
			'DAZ1', 'DAZ2', 'DAZ3', 'DAZ4', 'DDX3Y', 'EIF1AY', 'HSFY1', 'HSFY2',
			'KDM5D', 'NLGN4Y', 'PCDH11Y', 'PRKY',
			'PRY', 'PRY2', 'RBMY1A1', 'RBMY1B', 'RBMY1D', 'RBMY1E', 'RBMY1F',
			'RBMY1J', 'RPS4Y1', 'RPS4Y2', 'SRY', 'TBL1Y', 'TGIF2LY', 'TMSB4Y',
			'TSPY1', 'TSPY10', 'TSPY2', 'TSPY3', 'TSPY4', 'TSPY8', 'TSPY9P',
			'TXLNGY', 'USP9Y', 'UTY', 'VCY', 'VCY1B', 'XKRY', 'XKRY2', 'ZFY']

# maps gene names to multi-copy gene family name
XYMCFAMS = {
    'HSFX1':'HSFX', 'HSFX2':'HSFX', 
    'HSFY1':'HSFY', 'HSFY2':'HSFY',
    'TSPY1':'TSPY', 'TSPY2':'TSPY', 'TSPY3':'TSPY',
    'TSPY4':'TSPY', 'TSPY8':'TSPY', 'TSPY10':'TSPY', 'TSPY9P':'TSPY',
    'RBMY1A1':'RBMY', 'RBMY1B':'RBMY', 'RBMY1D':'RBMY',
    'RBMY1E':'RBMY', 'RBMY1F':'RBMY', 'RBMY1J':'RBMY',
    'VCX':'VCX', 'VCX2':'VCX', 'VCX3A':'VCX', 'VCX3B':'VCX',
    'VCY':'VCY', 'VCY1B':'VCY',
    'PRY':'PRY', 'PRY2':'PRY',
    'BPY2':'BPY2', 'BPY2B':'BPY2', 'BPY2C':'BPY2',
    'CDY1':'CDY', 'CDY1B':'CDY', 'CDY2A':'CDY', 'CDY2B':'CDY',
    'DAZ1':'DAZ', 'DAZ2':'DAZ', 'DAZ3':'DAZ', 'DAZ4':'DAZ',
    'XKRY':'XKRY', 'XKRY2':'XKRY'
}
XYMCFAM_NAMES = ['HSFY', 'TSPY', 'RBMY', 'VCY', 'PRY', 
				 'BPY2', 'CDY', 'DAZ', 'XKRY']

# X-Y pairs
XYPAIRS = [('AMELX', 'AMELY'), ('DDX3X', 'DDX3Y'), ('EIF1AX', 'EIF1AY'), 
           ('HSFX', 'HSFY'), ('KDM5C', 'KDM5D'), ('NLGN4X', 'NLGN4Y'),
           ('PCDH11X', 'PCDH11Y'), ('PRKX', 'PRKY'), ('RBMX', 'RBMY'), 
           ('RPS4X', 'RPS4Y1'), ('SOX3', 'SRY'), ('TBL1X', 'TBL1Y'),
           ('TGIF2LX', 'TGIF2LY'), ('TMSB4X', 'TMSB4Y'),
           ('TSPYL2', 'TSPY'), ('TXLNG', 'TXLNGY'), ('USP9X', 'USP9Y'),
           ('KDM6A', 'UTY'), ('VCX', 'VCY'), ('ZFX', 'ZFY')]
XYPAIRS_X = []
XYPAIRS_Y = []
for gx, gy in XYPAIRS:
    XYPAIRS_X.append(gx)
    XYPAIRS_Y.append(gy)

# the core group of the 9 most widely ancestral expressed X-Y pairs
XYPAIRS_X9 = ['DDX3X', 'EIF1AX', 'KDM5C', 'NLGN4X', 'PRKX', 'RPS4X',
              'USP9X', 'KDM6A', 'ZFX']
XYPAIRS_Y9 = ['DDX3Y', 'EIF1AY', 'KDM5D', 'NLGN4Y', 'PRKY', 'RPS4Y1',
              'USP9Y', 'UTY', 'ZFY']

# ...the above + TMSB4X/TMSB4Y and TBL1X/TBL1Y
XYPAIRS_X11 = ['DDX3X', 'EIF1AX', 'KDM5C', 'NLGN4X', 'PRKX', 'RPS4X',
               'TBL1X', 'TMSB4X', 'USP9X', 'KDM6A', 'ZFX']
XYPAIRS_Y11 = ['DDX3Y', 'EIF1AY', 'KDM5D', 'NLGN4Y', 'PRKY', 'RPS4Y1',
               'TBL1Y', 'TMSB4Y', 'USP9Y', 'UTY', 'ZFY']

# ...the above + PCDH11X/PCDH11Y
XYPAIRS_X12 = ['DDX3X', 'EIF1AX', 'KDM5C', 'NLGN4X', 'PCDH11X', 'PRKX', 'RPS4X',
               'TBL1X', 'TMSB4X', 'USP9X', 'KDM6A', 'ZFX']
XYPAIRS_Y12 = ['DDX3Y', 'EIF1AY', 'KDM5D', 'NLGN4Y', 'PCDH11Y', 'PRKY', 'RPS4Y1',
               'TBL1Y', 'TMSB4Y', 'USP9Y', 'UTY', 'ZFY']

# ...the above + TXLNG/TXLNGY
XYPAIRS_X13 = ['DDX3X', 'EIF1AX', 'KDM5C', 'NLGN4X', 'PCDH11X', 'PRKX', 'RPS4X',
               'TBL1X', 'TMSB4X', 'TXLNG', 'USP9X', 'KDM6A', 'ZFX']
XYPAIRS_Y13 = ['DDX3Y', 'EIF1AY', 'KDM5D', 'NLGN4Y', 'PCDH11Y', 'PRKY', 'RPS4Y1',
               'TBL1Y', 'TMSB4Y', 'TXLNGY', 'USP9Y', 'UTY', 'ZFY']


def sum_mc_families(data):
    """ Sum values of genes in multicopy gene families.
    
    Values of genes not in multicopy families will not be altered.
    
    Parameters
    ----------
    data : pandas.DataFrame
      genes x samples dataframe with gene names in index
    """
    allfams = {}
    for g in data.index:
        allfams[g] = XYMCFAMS[g] if g in XYMCFAMS else g
    data = data.groupby(by=lambda g: allfams[g], axis=0).sum()
    return data