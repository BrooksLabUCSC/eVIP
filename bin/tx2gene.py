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
def main( outDir=None, inDir=None, sampleList=None):

    '''
    main
    '''
    R.assign('inDir',inDir)
    R.assign('outdir',outDir)
    R.assign('sampleList',sampleList)


    # import
    from rpy2.robjects.packages import importr
    #kallisto processing libraries
    tximportData   = importr('tximportData')
    tximport   = importr('tximport')
    ensembldb   = importr('ensembldb')
    EnsDb_Hsapiens_v86   = importr('EnsDb.Hsapiens.v86')
    #deseq
    methods   = importr('methods')

    #transcripts to gene, used in tximport
    R('edb <- EnsDb.Hsapiens.v86')
    R('tx2gene = transcripts(edb , columns=c("tx_id", "gene_name"),return.type="DataFrame")')

    #locate kallisto files
    R('files <- file.path(inDir, sampleList, "abundance.tsv")')
    R('all(file.exists(files))')
    R('names(files)<- sampleList')

    #tximport conversion to gene
    R('txi.kallisto <- tximport(files, type = "kallisto",tx2gene = tx2gene, txOut = FALSE,ignoreTxVersion=TRUE, countsFromAbundance = "scaledTPM")')

    # R('print(head(txi.kallisto$counts))')
    # R('print(head(txi.kallisto$abundance))')

    #abundance is the combined gene count of the transcripts - $counts is...?
    R('counts = txi.kallisto$abundance')

    R('d <- cbind(rownames(counts), data.frame(counts, row.names=NULL))')
    R('colnames(d)[1] <- "#gene_id"')

    R('write.table(as.data.frame(d),file=outdir,row.names=FALSE,sep="\t",quote=FALSE)')


if __name__ == "__main__":
    main()
