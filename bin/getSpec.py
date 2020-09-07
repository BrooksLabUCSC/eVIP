########################################################################
# author: Alexis M. Thornton
#
########################################################################
import pandas as pd
import numpy as np
import itertools
from itertools import combinations
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.lines as mlines
import math
import os


def main(wt,mut,control,file_wt,file_mut,outdir):


    if not os.path.exists(outdir):
        os.makedirs(outdir)

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

    #to plot pval vs pval for each comparison
    if control_vs_wt_results_deseq2_counts.shape[0]>0 and control_vs_mutant_gene_level_sig.shape[0]>0:
        pvalplot(control_vs_wt_results_deseq2_counts,control_vs_mutant_gene_level_sig,mutspec,wtspec,common,mut,outdir)

    pvalhisto(control_vs_wt_gene_level_sig,control_vs_mutant_results_deseq2_counts,mut,outdir)
    FChisto(control_vs_wt_results_deseq2_counts,control_vs_mutant_results_deseq2_counts,mut,outdir)


    return mutspec,wtspec

def FChisto(control_vs_wt_results_deseq2_counts,control_vs_mutant_results_deseq2_counts,mut,outdir):

    plt.clf()

    plt.figure(figsize=(10,8))
    plt.hist(control_vs_wt_results_deseq2_counts["log2FoldChange"].dropna(),bins=70,alpha = 0.5,color='black',label="control vs WT")
    plt.hist(control_vs_mutant_results_deseq2_counts["log2FoldChange"].dropna(),bins=70,alpha = 0.5,color='red',label="control vs mutant")
    plt.legend()

    plt.title("Log2FC of DEG "+ mut)
    plt.xlabel('log2FC')
    plt.ylabel('count')
    plt.savefig(outdir+"/"+mut+"_logFC_histo.png", dpi=600, bbox_inches='tight')
    plt.clf()


def pvalhisto(control_vs_wt_results_deseq2_counts,control_vs_mutant_results_deseq2_counts,mut,outdir):

    plt.clf()

    plt.figure(figsize=(10,8))
    plt.hist(control_vs_wt_results_deseq2_counts["padj"].dropna(),bins=70,alpha = 0.5,color='black',label="control vs WT")
    plt.hist(control_vs_mutant_results_deseq2_counts["padj"].dropna(),bins=70,alpha = 0.5,color='red',label="control vs mutant")
    plt.legend()

    plt.title("Adjusted p-values "+ mut)
    plt.xlabel('adj pval')
    plt.ylabel('count')
    plt.savefig(outdir+"/"+mut+"_pval_histo.png", dpi=600, bbox_inches='tight')
    plt.clf()

def pvalplot(control_vs_wt_results_deseq2_counts,control_vs_mutant_results_deseq2_counts,mutspec,wtspec,common,mut,outdir):
        plot_df = pd.merge(control_vs_wt_results_deseq2_counts.rename(columns={"padj": "padj-WT","log2FoldChange":"log2FC-WT"})[["padj-WT","log2FC-WT"]],control_vs_mutant_results_deseq2_counts.rename(columns={"padj": "padj-mut","log2FoldChange":"log2FC-mut"})[["padj-mut","log2FC-mut"]] ,  left_index=True, right_index=True)

        plot_df["color"] = None
        plot_df['color'].loc[plot_df.index.isin(mutspec)] = '#C14242'
        plot_df['color'].loc[plot_df.index.isin(wtspec)] = '#1f77b4'
        plot_df['color'].loc[plot_df.index.isin(common)] =  '#000000'
        plot_df['color'] = plot_df['color'].fillna('#D3D3D3')


        plot_df["status"] = None
        plot_df['status'].loc[plot_df.index.isin(mutspec)] = 'mutation-specific'
        plot_df['status'].loc[plot_df.index.isin(wtspec)] = 'WT-specific'
        plot_df['status'].loc[plot_df.index.isin(common)] = 'DEG in both'
        plot_df['status'].loc[plot_df['color'] =='#D3D3D3'] = 'Not a DEG'

        smallest_num = np.nextafter(0, 1)

        plot_df["padj-WT"] = plot_df["padj-WT"]+smallest_num
        plot_df["padj-WT"] = plot_df["padj-WT"].to_frame().applymap(math.log10)*-1
        plot_df["padj-mut"] = plot_df["padj-mut"]+smallest_num
        plot_df["padj-mut"] = plot_df["padj-mut"].to_frame().applymap(math.log10)*-1

        groups = plot_df.groupby("status")

        plt.clf()
        plt.figure(figsize=(10,8))
        for name, group in groups:
            # print(group)

            plt.scatter(group["padj-WT"],group["padj-mut"], s=1,alpha=.2 , c=group["color"])

        #making legend
        mutspecific = mlines.Line2D([], [], color='#C14242', marker='o', linestyle='None',markersize=3, label='mutation-specific DEG')
        wtspecific = mlines.Line2D([], [], color='#1f77b4', marker='o', linestyle='None', markersize=3, label='WT-specific DEG')
        common = mlines.Line2D([], [], color='#000000', marker='o', linestyle='None', markersize=3, label='DEG in both')
        notDEG = mlines.Line2D([], [], color='#D3D3D3', marker='o', linestyle='None',markersize=3, label='Not a DEG')

        plt.legend(handles=[wtspecific, mutspecific, common, notDEG],frameon=False,bbox_to_anchor=(1.005, 1), loc=2, borderaxespad=0.)

        plt.title("Adjusted p-values distribution "+ mut)
        plt.xlabel('control vs WT : -log10(adj pval)')
        plt.ylabel('control vs mutant : -log10(adj pval)')
        plt.savefig(outdir+"/"+mut+"_pval_plot.png", dpi=600, bbox_inches='tight')
        plt.clf()

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
