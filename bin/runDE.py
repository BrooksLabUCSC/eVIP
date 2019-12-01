########################################################################
# Adapted from runDE.py in FLAIR
# https://github.com/BrooksLabUCSC/flair
# Original author: Cameron M. Soulette
#
########################################################################

import os, sys
import pandas as pd
import numpy as np

#supressing rpy2 warnings
import warnings
from rpy2.rinterface import RRuntimeWarning
warnings.filterwarnings("ignore", category=RRuntimeWarning)

from rpy2 import robjects
from rpy2.robjects import r,pandas2ri, Formula
from rpy2.robjects.lib import grid
pandas2ri.activate()
R = robjects.r

# main
def main(group1=None, group2=None, outDir=None, inDir=None, formula=None):

    '''
    main
    '''

    R.assign('inDir',inDir)
    R.assign('outdir',outDir)

    R.assign('group1',group1)
    R.assign('group2',group2)

    print("Running DeSeq2....")
    print(group1 +" vs "+ group2)

    # import
    from rpy2.robjects.packages import importr
    #kallisto processing libraries
    tximportData   = importr('tximportData')
    tximport   = importr('tximport')
    ensembldb   = importr('ensembldb')
    EnsDb_Hsapiens_v86   = importr('EnsDb.Hsapiens.v86')
    #deseq
    methods   = importr('methods')
    deseq     = importr('DESeq2')

    #transcripts to gene, used in tximport
    R('edb <- EnsDb.Hsapiens.v86')
    R('tx2gene = transcripts(edb , columns=c("tx_id", "gene_name"),return.type="DataFrame")')

    # import formula
    formulaDF     = pd.read_csv(formula,header=0, sep="\t")

    samples =  formulaDF.samples.tolist()
    R.assign('samples',samples)

    sampleTable = pandas2ri.py2ri(formulaDF)
    R.assign('sampleTable',sampleTable)

    #locate kallisto files
    R('files <- file.path(inDir, samples, "abundance.h5")')
    R('all(file.exists(files))')

    #tximport conversion to gene
    R('txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene, txOut = FALSE,ignoreTxVersion=TRUE)')
    R('rownames(sampleTable) <- samples')

    #DESeq
    R('dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, ~condition)')

    # R('colData(dds)$condition<-factor(colData(dds)$condition, levels=c(group1,group2))')

    R('dds_<-DESeq(dds)')
    R('res<-results(dds_)')
    R('res<-res[order(res$padj),]')

    # writing deseq2 results to a file
    Out = os.path.join(outDir, "%s_v_%s_deseq2_results.tsv"  % (group1,group2))
    R.assign('Out',Out)

    R('write.csv(as.data.frame(res),file=Out)')


if __name__ == "__main__":
    main()
