########################################################################
# author: Alexis M. Thornton
#
########################################################################
import pandas as pd
import numpy as np
from itertools import combinations


def main(wt,mut,control,file_wt,file_mut):

    print "\n"
    print "Getting gene lists for",mut,"..."

    control_vs_wt_results_deseq2_counts = pd.read_csv(file_wt, header = 0, index_col=0)
    control_vs_mutant_results_deseq2_counts = pd.read_csv(file_mut, header = 0, index_col=0)

    control_vs_wt_gene_level_sig = control_vs_wt_results_deseq2_counts[(control_vs_wt_results_deseq2_counts['padj']<=0.05) & (abs(control_vs_wt_results_deseq2_counts['log2FoldChange'])>0.5)]
    control_vs_wt_gene_level_sig_genes = control_vs_wt_gene_level_sig.index.tolist()

    control_vs_mutant_gene_level_sig = control_vs_mutant_results_deseq2_counts[(control_vs_mutant_results_deseq2_counts['padj']<=0.05)& (abs(control_vs_mutant_results_deseq2_counts['log2FoldChange'])>0.5)]
    control_vs_mutant_gene_level_sig_genes = control_vs_mutant_gene_level_sig.index.tolist()


    mutspec, wtspec, common = getMutSpec(control_vs_mutant_gene_level_sig_genes, control_vs_wt_gene_level_sig_genes)

    print "mutation specific genes:", len(mutspec)
    print "WT specific genes:", len(wtspec)
    print "common DEG genes:", len(common)

    return mutspec,wtspec

def getMutSpec(control_vs_mut_genes, control_vs_WT_genes):

    #function to get mutation specific genes, WT specific genes,
    #and commonly differentially expressed genes

    """
    inputs = different list of significant genes
    outputs = differentially expressed genes that are
    specific to the mutation, WT, or commonly differentially expressed
    """

    data = dict(
    controlvsmut = set(control_vs_mut_genes),
    controlvsWT = set(control_vs_WT_genes),
    )

    variations = {}
    for i in range(len(data)):
        for v in combinations(data.keys(),i+1):
            vsets = [ data[x] for x in v ]
            variations[tuple(sorted(v))] = reduce(lambda x,y: x.intersection(y), vsets)

    #Mut only
    mutspec = list(set(variations[ ('controlvsmut',)]) - set(variations[ ('controlvsWT',)]))

    #genes in middle
    common = set(control_vs_WT_genes).intersection(set(control_vs_mut_genes))

    # WT only
    wtspec= list(set(variations[ ('controlvsWT',)]) - set(variations[ ('controlvsmut',)]))

    return mutspec, wtspec, common


if __name__ == "__main__":
    main()
