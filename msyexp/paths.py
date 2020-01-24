""" paths to input files """
import os

##### ~~~ EDIT HERE ~~~ #######################################################

# local directory containing all input files downloaded from Zenodo
# (https://doi.org/10.5281/zenodo.3627110)
DATADIR =    # your path here

# local directory where you would like all output files to be written;
# each notebook's files will be saved within a notebook-specific subdirectory
NBOUTDIR =   # your path here

##### ~~~ EDIT HERE ~~~ #######################################################


##### ~~~ DO NOT CHANGE BELOW THIS POINT ~~~ ##################################

if not os.path.exists(NBOUTDIR):
	os.mkdir(NBOUTDIR)

# GTEx samples --> kallisto
EXP_KSTO = '{}/gtex.ksto.gen24c.prefilt_pass1.postfilt_pass.pcgenes.tpm.h5'.format(DATADIR)
# adjusted expression levels (RIN, PMI, intronic read mapping) for male samples only
EXP_KSTO_ADJ_M = '{}/gtex.ksto.gen24c.prefilt_pass1.postfilt_pass.pcgenes.males.adj_tpm.h5'.format(DATADIR)
# adjusted expression levels, both males and females
EXP_KSTO_ADJ_MF = '{}/gtex.ksto.gen24c.prefilt_pass1.postfilt_pass.pcgenes.both.adj_tpm.h5'.format(DATADIR)

# unadjusted transcript TPMs
EXP_TX_KSTO = '{}/gtex.ksto.gen24c.prefilt_pass1.postfilt_pass.transcripts.tpm.h5'.format(DATADIR)

# GTEx samples --> RNA-SeQC (GTEx Consortium analysis)
EXP_RNASEQC = '{}/gtex.rnaseqc.gen19.pcgenes.tpm.h5'.format(DATADIR)

## HPA samples --> kallisto
EXP_HPA_KSTO = '{}/hpa.ksto.gen24c.pcgenes.tpm.txt'.format(DATADIR)

# sample/donor metadata
METADATA_GTEX = '{}/gtex_metadata_subset.txt'.format(DATADIR)
METADATA_HPA = '{}/hpa.metadata.with_gtex_tissues.txt'.format(DATADIR)

# gene/transcript annotation data (based on GENCODE v24 GTF)
ANNOFILE = '{}/gencode.v24.annotation.basic_ccds_nopar.gene_tx_annotable.txt'.format(DATADIR)

# results from RNA-seq simulations
SIM_Y0 = '{}/sim_y0.tpm.txt'.format(DATADIR)
SIM_Y1_X2 = '{}/sim_y1_x2.tpm.txt'.format(DATADIR)
SIM_Y5_X10 = '{}/sim_y5_x10.tpm.txt'.format(DATADIR)
SIM_Y5_XUNCORR = '{}/sim_y5_xuncorr.tpm.txt'.format(DATADIR)

# mappability data
MAPP_DATA_V19 = '{}/gen19.pcgenes_longesttx.76bp_tiled_reads.mappability.txt'.format(DATADIR)

# results from limma sex bias analysis
LIMMA_SEXBIAS = '{}/xypair_transcript_sexbias.limma_xycombined.txt'.format(DATADIR)

# Ludwig et al. 2016 miRNA expression data
MIR1_EXP = '{}/ludwig2016_mir1_expression.txt'.format(DATADIR)

# mass spec data
MASS_SPEC_PEPTIDES = '{}/mass_spec.peptide_quant.h5'.format(DATADIR)
MASS_SPEC_PROTEINS = '{}/mass_spec.protein_quant.h5'.format(DATADIR)
MASS_SPEC_METADATA = '{}/mass_spec_metadata.txt'.format(DATADIR)

# Brawand/Merkin data
BRAMERK_META = "{}/brawand_merkin.metadata.txt".format(DATADIR)
BRAMERK_CHICKEN = '{}/brawand_merkin.chicken.transcripts.tpm.txt'.format(DATADIR)
BRAMERK_CHIMP = '{}/brawand.chimp.transcripts.tpm.txt'.format(DATADIR)
BRAMERK_MOUSE = '{}/brawand_merkin.mouse.transcripts.tpm.txt'.format(DATADIR)
BRAMERK_RHESUS = '{}/brawand_merkin.rhesus.transcripts.tpm.txt'.format(DATADIR)

# amino-acid sequences for X-Y pairs
UNIPROT_XYSEQS = '{}/uniprot_xy_protein_seqs.fa'.format(DATADIR)


