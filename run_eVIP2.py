# Author: Alexis M Thornton

#!/usr/bin/python
import sys
import argparse
import os
import errno
import csv
import itertools
import rpy2.robjects as robjects
import json
import pandas as pd
from pyGeno.tools.parsers.GTFTools import GTFFile

#importing eVIP
from bin import eVIP_corr
from bin import eVIP_predict
from bin import eVIP_sparkler
from bin import eVIP_viz
from bin import eVIP_compare
from bin import eVIPP_sparkler

#eVIP2
from bin import filterGeneExpressionTable
from bin import runDE
from bin import getSpec
from bin import eVIPPspec
from bin import combine_sparklers
from bin import upset_plot
from bin import tx2gene

########
# MAIN #
########
def main(infile=None, zscore_gct = None, out_directory=None, sig_info =None,
        c=None, r=None, num_reps=None,ie_filter=None,ie_col=None, i=None,
        allele_col=None, conn_null=None, conn_thresh=None,
        mut_wt_rep_rank_diff=None, use_c_pval=None, cell_id=None, plate_id=None,
        ref_allele_mode=None,x_thresh=None, y_thresh=None, annotate=None,
        by_gene_color=None, pdf=None, xmin=None,xmax=None, ymin=None, ymax=None,
        viz_ymin=None, viz_ymax=None, corr_val=None):

    parser = argparse.ArgumentParser()

    #from filter gene expression table
    parser.add_argument("--min_tpm", help = """ When filtering the gene expression 
                table, this value is the minimum TPM value for a given
                gene. If the gene is expressed below this level in all samples,
                the gene is removed from the table. DEFAULT=1""",
                default=1,type=float)

    #from corr
    parser.add_argument("--infile", help="""Input txt file (filtered and
                log transformed data).""")
    parser.add_argument("--input_gene_tpm",
                help="Gene tpm table input for eVIP overall prediction")
    parser.add_argument("-zscore_gct", help="""Zscore input gct file (use
                instead of --infile)""")
    parser.add_argument("--out_directory",required=True, help="""Path to directory
                for eVIP output files""")
    #from compare
    parser.add_argument("--sig_info",required=True, help = """sig info file with
                gene information and distil information""")
    parser.add_argument("-c",required=True, help = """.grp file containing
                allele names of control perturbations. If this file is given,
                a null will be calculated from these""")
    parser.add_argument("-r", required=True, help = """File explicitly
                indicating which comparisons to do. Assumes the file has a
                header and it is ignored. The first column is the reference
                allele and second column is test allele. If this file is not
                given, then the reference alleles are assumed to be WT and
                inferred from the allele names.""")
    parser.add_argument("--num_reps",required=True, help = """Number of
                replicates expected for each allele.""")
    parser.add_argument("--ie_filter", help = """Threshold for infection
                efficiency in L1000. Any wildtype or mutant alleles having an ie below
                this threshold, will be removed""")
    parser.add_argument("--ie_col", help = """Name of the column in the sig_info
                file with infection efficiency information.""")
    # parser.add_argument("-i", help = "Number of iterations to run. DEF=1000")
    parser.add_argument("--allele_col", default = "allele", help = """Column name
                in sig_info file that indicates the allele names.DEF=allele""")
    parser.add_argument("--conn_null", help = """ Optional file containing
                connectivity null values from a previous run. Should end
                in _conn_null.txt""")
    #from predict
    parser.add_argument("--conn_thresh",help = """P-value threshold for
                connectivity vs null. DEFAULT=0.1""",
                default=0.1,type=float)
    parser.add_argument("--mut_wt_rep_thresh",
                help = """P-value threshold for comparison of WT and mut
                robustness. DEFAULT=0.1""",
                default=0.1, type=float)
    parser.add_argument("--disting_thresh", help = """P-value threshold that
                tests if mut and wt reps are indistinguishable from each other.
                DEFAULT=0.1""",
                default=0.1,type=float)
    parser.add_argument("--mut_wt_rep_rank_diff", help = """The minimum
                difference in median rankpoint WT and mut to consider a
                difference. DEF=0""", default=0, type=float)
    parser.add_argument("--use_c_pval", action ="store_true",
                help = "Will use corrected p-value instead of raw p-val")
    parser.add_argument("--cell_id",
                help = """Optional: Will only look at signatures from this cell
                line. Helps to filter sig_info file.""")
    parser.add_argument("--plate_id",
                help = "Optional: Will only look at signatures from this plate")
    parser.add_argument("--cond_max_diff_thresh",
                help = """Threshold for maximum difference between condition
                correlation medians when determining if variant is not neutral.
                Default = 0.2""",
                type=float,default=0.2)

    #from sparkler
    parser.add_argument("--ref_allele_mode", action ="store_true",
                help = """Sparkler+Viz: Instead of organizing plots by gene,
                will use the wt column to determine what are the reference
                alleles.""" )
    parser.add_argument("--x_thresh" ,
                help = "Sparkler: Threshold of significance",
                default=1.3,type=float)
    parser.add_argument("--y_thresh",
                help = "Sparkler: Threshold of impact direction",
                default=1.3,type=float)
    parser.add_argument("--annotate", action ="store_true",
                help = "Sparkler: Will add allele labels to points.")
    parser.add_argument("--by_gene_color",
                help = """Sparkler: File containing labels and colors for
                gene-centric plot.""")
    parser.add_argument("--pdf",
                help = """Sparkler + Viz: Will print plots in pdf format instead
                of png.""")
    parser.add_argument("--xmin",
                help = "Sparkler: Min value of x-axis. DEF=0",
                type=float,default=0)
    parser.add_argument("--xmax",
                help = "Sparkler: Max value of x-axis. DEF=4",
                type=float,default=4)
    parser.add_argument("--ymin",
                help = "Sparkler: Min value of y-axis. DEF=-3",
                type=float,default=-3)
    parser.add_argument("--ymax",
                help = "Sparkler: Min value of y-axis. DEF=3",
                type=float,default=3)
    #from viz
    parser.add_argument("--viz_ymin",
                help = "Viz: Minimum y-value of rep value. DEF=-1",
                type=float,default=-1)
    parser.add_argument("--viz_ymax",
                help = "Viz: Maximum y-value of rep value. DEF=1",
                type=float,default=1)
    parser.add_argument("--corr_val",
                help = """Viz: String used to label the correlation value.
                DEF= 'spearman' """, default = "spearman")
    #eVIPP
    parser.add_argument("--eVIPP", action ="store_true",
                help="""Use this option when doing pathway analysis, must also
                have gmt or JSON file """)
    parser.add_argument("--JSON",
                help= """JSON file created by create_pathway_JSON.py. Contains
                dictionary of pathways and the associated ids""")
    parser.add_argument("--gmt", help= "Gene set file in .gmt format")
    parser.add_argument("--min_genes",
                help = """Minimum amount of pathway genes found in data to run
                eVIPP on. DEF = 10""",default=10, type=float)
    parser.add_argument("--viz_off", action ="store_true",
                help = "Will not perform eVIP viz step")
    parser.add_argument("-sparkler_off", action ="store_true",
                help = "Will not perform eVIP sparkler step")

    #run_eVIP2
    parser.add_argument("--input_dir",
                help="Path to directory of kallisto outputs")
    parser.add_argument("--gtf",
                help="Gtf file used to convert transcript counts to gene counts")
    parser.add_argument("--control",
                required=False,
                help="""If multiple controls in the controls file, designate
                which to use for deseq2""")
    parser.add_argument("--tx2gene",
                action ="store_true",required=False,
                help="""Use tximport for transcript to gene conversion when
                using -input_dir""")

    global args
    args = parser.parse_args()

    #make eVIP output directory
    global out_dir
    out_dir = args.out_directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    #must have atleast one kind of input
    input_types =[args.input_dir , args.input_gene_tpm, args.infile, args.zscore_gct]
    if all(v is None for v in input_types):
        print("Error:Input data is missing")
        sys.exit()
    #If runnig eVIPP
    if args.eVIPP :
        if all(v is None for v in [args.JSON,args.gmt]):
            print("Error: Must input gmt or JSON pathway file when running eVIPP")
            sys.exit()


    #############################################################################
    ### running overall eVIP from kallisto outputs

    if args.input_dir:

        #combining kallisto abundance files into one file
        combined_kallisto_transcript_df,all_samples = kallisto_process()

        if args.tx2gene:
            # new version using tximport to combine kallisto into gene counts
            # using kallisto directories
            tx2gene.main(outDir=args.out_directory+"/kallisto_files/combined_kallisto_abundance_genes.tsv",
                        inDir=args.input_dir, sampleList = all_samples )

        else:
            #combining to gene level the original way
            transcript_to_gene_counts(combined_kallisto_transcript_df)

        #filtering out low expressed genes and doing log2 transformation
        print("Filtering out low expressed genes and doing log2 transformation...")
        filterGeneExpressionTable.main(in_table=args.out_directory+"/kallisto_files/combined_kallisto_abundance_genes.tsv",
                        out_table=args.out_directory+"/kallisto_files/combined_kallisto_abundance_genes_filtered_transformed.tsv",
                        x = 1,l=True,reformat_gene = None,fpkms = None,
                        min_val = args.min_tpm, min_fold_fpkm = None)

        eVIP_infile_path = args.out_directory+"/kallisto_files/combined_kallisto_abundance_genes_filtered_transformed.tsv"

    #if input gene tpm
    if args.input_gene_tpm:
        
        filterGeneExpressionTable.main(in_table=args.input_gene_tpm,
                    out_table=args.out_directory+"/combined_kallisto_abundance_genes_filtered_transformed.tsv",
                    x = 1,l=True,reformat_gene = None,fpkms = None,
                    min_val = args.min_tpm, min_fold_fpkm = None)

        eVIP_infile_path = args.out_directory+"/combined_kallisto_abundance_genes_filtered_transformed.tsv"


    if args.infile: 
        eVIP_infile_path = args.infile

    #run eVIP overall
    overall_eVIP_dir = args.out_directory + "/eVIP_out"

    if not os.path.exists(overall_eVIP_dir):
        os.makedirs(overall_eVIP_dir)


    print("Running eVIP for overall function...")
    run_eVIP(eVIP_infile_path,
            None, overall_eVIP_dir, args.sig_info, args.c, args.r, args.num_reps,
            args.ie_filter, args.ie_col, None, args.allele_col, args.conn_null,
            args.conn_thresh,args.mut_wt_rep_rank_diff, args.use_c_pval,
            args.cell_id, args.plate_id, args.ref_allele_mode,args.x_thresh,
            args.y_thresh, args.annotate, args.by_gene_color, args.pdf, args.xmin,
            args.xmax, args.ymin, args.ymax, args.viz_ymin, args.viz_ymax,
            args.corr_val)


    #############################################################################
    ### eVIP Pathways

    if args.eVIPP and args.input_dir:

        print ("Running eVIP Pathways...")

        if not os.path.exists(args.out_directory + "/eVIPP_out"):
            os.makedirs(args.out_directory + "/eVIPP_out")

        #make dictionary of replicates to condition from sig info
        cond_to_rep_dict = condition_to_replicates(args.sig_info)

        #for each mutant, do control vs mutant comparison

        #finding which control to  use
        controls_list = []
        with open(args.c) as control_file:
            controls_list = [line.strip() for line in control_file]

        #if more than one use the -control argument
        if len(controls_list) > 1:
            deseq_control = args.control
        else:
            deseq_control = controls_list[0]

        if not deseq_control:
            print("Error: Need deseq control. Set with --control ")
            sys.exit()

        comparisons_df = pd.read_csv(args.r, delim_whitespace=True)
        comparisons_dict = comparisons_df.to_dict('records')

        deseq2_muts(comparisons_dict,cond_to_rep_dict,deseq_control)

        #find all WTs in the comparison dict
        wts = list(set([i['wt'] for i in comparisons_dict]))

        #run deseq on all wts
        deseq2_wts(comparisons_dict,cond_to_rep_dict,deseq_control,wts)

        #for each mutation
        for i in comparisons_dict:
            wt = i['wt']
            mut = i['mutant']

            file_mut = args.out_directory+"/deseq2/"+deseq_control+"_vs_"+mut+"/"+deseq_control+"_v_"+mut+"_deseq2_results.csv"
            file_wt = args.out_directory+"/deseq2/"+deseq_control+"_vs_"+wt+"/"+deseq_control+"_v_"+wt+"_deseq2_results.csv"

            #get mutation specific and wt specific genes
            mutspec,wtspec = getSpec.main(wt,mut,deseq_control,file_wt,
                                file_mut,args.out_directory+"/deseq2/figures")

            eVIP_gene_expression = pd.read_csv(args.out_directory+"/kallisto_files/combined_kallisto_abundance_genes_filtered_transformed.tsv", sep = "\t")

            ######################################################################
            #making eVIP files
            eVIPP_files = args.out_directory + "/eVIPP_out/"+mut+"/eVIP_files"

            if not os.path.exists(eVIPP_files):
                os.makedirs(eVIPP_files)

            #making new sig_info
            with open(args.sig_info) as sig_info, open(eVIPP_files+"/"+mut+"_sig.info","w+") as spec_sig_info:
                for line in sig_info:
                    if line.split()[1] in [wt,mut,deseq_control,"sig_id"]:
                        spec_sig_info.write(line)

            #making new comparisons file
            with open(args.r) as r, open(eVIPP_files+"/"+mut+"_comparisons.tsv","w+") as spec_comparisons:
                spec_comparisons.write(("\t").join(["wt","mutant"]))
                spec_comparisons.write("\n")

                for line in r:
                    if line.split()[1] == mut:
                        spec_comparisons.write(line)

            #################################################################
            #mutation specific

            #subset gene expression file - mutation specific
            eVIP_gene_expression_mutspec = eVIP_gene_expression[eVIP_gene_expression['#gene_id'].isin(mutspec)]

            #how many mutspec genes were in the file?
            print "Number of mutation-specific genes present in the filtered kallisto file:",eVIP_gene_expression_mutspec.shape[0], "of" ,len(mutspec)
            print ("\n")

            mutspec_infile = args.out_directory+"/kallisto_files/combined_kallisto_abundance_genes_filtered_transformed_"+mut+"_mutspec.tsv"
            eVIP_gene_expression_mutspec.to_csv(mutspec_infile,sep="\t",index=False)


            ######################################################################
            #run eVIPP on mutation specific

            eVIPP_mutspec_out = args.out_directory + "/eVIPP_out/"+mut+"/mutation_specific"

            if not os.path.exists(eVIPP_mutspec_out):
                os.makedirs(eVIPP_mutspec_out)


            eVIPPspec.main(eVIPP_mutspec_out,args.JSON,args.gmt,args.min_genes,
                        mutspec_infile,eVIPP_files+"/"+mut+"_sig.info",
                        args.c, eVIPP_files+"/"+mut+"_comparisons.tsv",
                        args.num_reps,args.ie_filter, args.ie_col, None,
                        args.allele_col, args.conn_null, args.conn_thresh,
                        args.mut_wt_rep_rank_diff, args.use_c_pval, args.cell_id,
                        args.plate_id, args.ref_allele_mode,args.x_thresh,
                        args.y_thresh, args.annotate, args.by_gene_color,
                        args.pdf, args.xmin,args.xmax, args.ymin, args.ymax,
                        args.viz_ymin, args.viz_ymax, args.corr_val,
                        args.mut_wt_rep_thresh,args.disting_thresh,
                        args.sparkler_off,args.viz_off,args.cond_max_diff_thresh)


            if os.path.exists(eVIPP_mutspec_out+"/eVIPP_combined_predict_files.txt"):
                upset_plot.run(args.JSON,args.gmt,mutspec_infile,
                        eVIPP_mutspec_out+"/eVIPP_combined_predict_files.txt",
                        eVIPP_mutspec_out+"/eVIPP_gene_overlap.png")


            ######################################################################
            # wt specific

            #subset gene expression file - wt specific
            eVIP_gene_expression_wtspec = eVIP_gene_expression[eVIP_gene_expression['#gene_id'].isin(wtspec)]

            #how many mutspec genes were in the file?
            print "Number of wt-specific genes present in the filtered kallisto file:",eVIP_gene_expression_wtspec.shape[0], "of" ,len(wtspec)
            print ("\n")

            wtspec_infile = args.out_directory+"/kallisto_files/combined_kallisto_abundance_genes_filtered_transformed_transformed_"+mut+"_wtspec.tsv"
            eVIP_gene_expression_wtspec.to_csv(wtspec_infile,sep="\t",index=False)

            ######################################################################
            #run eVIPP on wt specific

            eVIPP_wtspec_out = args.out_directory + "/eVIPP_out/"+mut+"/wt_specific"

            if not os.path.exists(eVIPP_wtspec_out):
                os.makedirs(eVIPP_wtspec_out)

            eVIPPspec.main(eVIPP_wtspec_out,args.JSON,args.gmt,args.min_genes,
                        wtspec_infile,eVIPP_files+"/"+mut+"_sig.info", args.c,
                        eVIPP_files+"/"+mut+"_comparisons.tsv", args.num_reps,
                        args.ie_filter, args.ie_col, None, args.allele_col,
                        args.conn_null, args.conn_thresh,args.mut_wt_rep_rank_diff,
                        args.use_c_pval, args.cell_id, args.plate_id,
                        args.ref_allele_mode,args.x_thresh, args.y_thresh,
                        args.annotate, args.by_gene_color, args.pdf, args.xmin,
                        args.xmax, args.ymin, args.ymax, args.viz_ymin,
                        args.viz_ymax, args.corr_val,args.mut_wt_rep_thresh,
                        args.disting_thresh,args.sparkler_off,args.viz_off,
                        args.cond_max_diff_thresh)

            if os.path.exists(eVIPP_wtspec_out+"/eVIPP_combined_predict_files.txt"):
                upset_plot.run(args.JSON,args.gmt,wtspec_infile,
                                eVIPP_wtspec_out+"/eVIPP_combined_predict_files.txt",
                                eVIPP_wtspec_out+"/eVIPP_gene_overlap.png")


        #####
        # combine all eVIPP mutation-specific and WT-specific sparklers into a single report
        combine_sparklers.run(args.out_directory + "/eVIPP_out",args.out_directory + "/eVIPP_out/all_eVIPP_sparklers.png")

        #####
        # use all mutation eVIPP predict files to make an output with each mutation and pathway prediction
        #list each directory(mutation)
        eVIPP_mut_paths = [args.out_directory + "/eVIPP_out/"+dI for dI in os.listdir(args.out_directory + "/eVIPP_out/") if os.path.isdir(os.path.join(args.out_directory + "/eVIPP_out/",dI))]

        #getting wt and mut specific eVIPP predict outputs and removing ones that dont exist
        mut_spec_files = [ i+"/mutation_specific/eVIPP_combined_predict_files.txt" for i in eVIPP_mut_paths if os.path.isfile(i+"/mutation_specific/eVIPP_combined_predict_files.txt")]
        wt_spec_files = [ i+"/wt_specific/eVIPP_combined_predict_files.txt" for i in eVIPP_mut_paths if os.path.isfile(i+"/wt_specific/eVIPP_combined_predict_files.txt")]

        if len(mut_spec_files)> 0:
            mut_spec_combined_df =  make_combined_pathway_df(mut_spec_files)
            mut_spec_combined_df.to_csv(args.out_directory + "/eVIPP_out/all_mutation_specific_eVIPP_summary.txt",sep="\t")

        if len(wt_spec_files)> 0:
            wt_spec_combined_df =  make_combined_pathway_df(wt_spec_files)
            wt_spec_combined_df.to_csv(args.out_directory + "/eVIPP_out/all_wt_specific_eVIPP_summary.txt",sep="\t")

    # Create the html interactive visualization with python3
    if args.use_c_pval :
        cmd = "python3 ./bin/generate_visualizations.py --inp_dir " +args.out_directory + " --use_c_pval "
        os.system(cmd)
    else:
        cmd = "python3 ./bin/generate_visualizations.py --inp_dir " +args.out_directory
        os.system(cmd)

#############
# FUNCTIONS #
#############

def make_combined_pathway_df(file_list):
    df = pd.concat([pd.read_csv(f,sep =  "\t") for f in file_list])
    df = df[["Pathway","mut","prediction"]]
    df = df.set_index("Pathway")
    df = df.pivot(columns='mut')
    df.columns = df.columns.droplevel()

    return df

def deseq2_wts(comparisons_dict,cond_to_rep_dict,deseq_control,wts):
    for wt in wts:
        comparison_dir = args.out_directory+"/deseq2/"+deseq_control+"_vs_"+wt
        if not os.path.exists(comparison_dir):
            os.makedirs(comparison_dir)

        #make formula df for running DeSeq2
        df = pd.DataFrame(columns=['samples','condition'])

        #adding control replciates
        for j in cond_to_rep_dict[deseq_control]:
            new_row = {'samples':j,'condition' :deseq_control}
            df = df.append(new_row, ignore_index=True)

        #adding wt replicates
        for k in cond_to_rep_dict[wt]:
            new_row = {'samples':k,'condition' :wt}
            df = df.append(new_row, ignore_index=True)

        #save df
        formula_file = comparison_dir +"/formula.tsv"
        df.to_csv(formula_file, sep="\t", index=False)

        #run DESeq 2
        runDE.main(group1=deseq_control, group2=wt, outDir=comparison_dir,
                    inDir=args.input_dir, formula=formula_file)

def deseq2_muts(comparisons_dict,cond_to_rep_dict,deseq_control):

    #for each mutation
    for i in comparisons_dict:
        # wt = i['wt']
        mut = i['mutant']

        comparison_dir = args.out_directory+"/deseq2/"+deseq_control+"_vs_"+mut
        if not os.path.exists(comparison_dir):
            os.makedirs(comparison_dir)

        #make formula df for running DeSeq2
        df = pd.DataFrame(columns=['samples','condition'])

        #adding control replciates
        for j in cond_to_rep_dict[deseq_control]:
            new_row = {'samples':j,'condition' :deseq_control}
            df = df.append(new_row, ignore_index=True)

        #adding mutant replicates
        for k in cond_to_rep_dict[mut]:
            new_row = {'samples':k,'condition' :mut}
            df = df.append(new_row, ignore_index=True)

        #save df
        formula_file = comparison_dir +"/formula.tsv"
        df.to_csv(formula_file, sep="\t", index=False)

        #run DESeq 2
        runDE.main(group1=deseq_control, group2=mut, outDir=comparison_dir,
                    inDir=args.input_dir, formula=formula_file)

def condition_to_replicates(sig_info_file):
    cond_to_rep = {}
    with open(args.sig_info) as sig_info:
        for line in sig_info:
            if line.startswith("distil_id"):
                continue
            reps = line.split("\t")[0].split("|")
            cond = line.split("\t")[1]
            cond_to_rep[cond]=reps
    return cond_to_rep

def transcript_to_gene_counts(transcript_df):
    print "Calculating gene counts from transcript counts..."

    #need to remove decimal from transcript id
    transcript_df_nodec = transcript_df
    transcript_df_nodec.index = transcript_df.index.str[:15]

    #making a dictionary of transcript id to gene name

    gtf = GTFFile(args.gtf)
    out = open(args.out_directory+"/kallisto_files/transcripts_genes.tsv", "w")
    transcript2gene = {}

    for line in gtf:
        try:
            transcript2gene[line['transcript_id']]=line["gene_name"]
            out.write(line['transcript_id'] +"\t"+ line["gene_name"]+ "\n")
        #when there is no transcript id in gtf line /no gene name for the transcript id
        except:
            pass

    #combining transcript to gene counts

    #make new df for gene counts
    gene_data = pd.DataFrame(columns=transcript_df_nodec.columns)
    gene_data.columns.name = '#gene_id'

    #for each column/sample
    for i in transcript_df_nodec:
        estimated_counts_per_gene={}
        #for each row in this column
        for transcript, row in transcript_df_nodec[i].iteritems():
            #count dict
            try:
                gene = transcript2gene[transcript]
                if gene in estimated_counts_per_gene:
                    estimated_counts_per_gene[gene] += row
                else:
                    estimated_counts_per_gene[gene] = row
            #pass if transcript not in conversion dictionary
            except:
                pass

        gene_data[i] = pd.Series(estimated_counts_per_gene)

    gene_data.to_csv(args.out_directory+"/kallisto_files/combined_kallisto_abundance_genes.tsv",
                            sep="\t", index_label="#gene_id")

def kallisto_process():
    print "Processing Kallisto files..."
    #get all sample names , directory names must match names in sig_info
    samples_sig = []
    with open(args.sig_info) as sig_info:
        for line in sig_info:
            for rep in line.split("\t")[0].split("|"):
                samples_sig.append(rep)
    samples = samples_sig[1:]

    #dict to match tsvs to sample abundance.tsv paths
    sample_to_tsv = {}
    for sample in samples:
        sample_to_tsv[sample] = args.input_dir+"/"+sample+"/abundance.tsv"

    #combining df into one
    data=pd.DataFrame()
    for i in sample_to_tsv.keys():
        table=pd.read_csv(sample_to_tsv[i], sep="\t", index_col=0)
        #keep only tpm column
        table = table[["tpm"]]
        #rename tpm column to sample name
        table = table.rename(columns={"tpm": i})
        #concat
        data=pd.concat([data,table], axis=1)

    if not os.path.exists(args.out_directory+"/kallisto_files"):
        os.makedirs(args.out_directory+"/kallisto_files")

    data.to_csv(args.out_directory+"/kallisto_files/combined_kallisto_abundance.tsv",
                        sep="\t")

    return data , samples

def run_eVIP(infile=None, zscore_gct = None, out_directory=None, sig_info =None,
            c=None, r=None, num_reps=None,ie_filter=None,ie_col=None, i=None,
            allele_col=None, conn_null=None, conn_thresh=None,
             mut_wt_rep_rank_diff=None, use_c_pval=None, cell_id=None,
             plate_id=None, ref_allele_mode=None,x_thresh=None, y_thresh=None,
             annotate=None, by_gene_color=None, pdf=None, xmin=None,xmax=None,
             ymin=None, ymax=None, viz_ymin=None, viz_ymax=None, corr_val=None):

    #different sig_gctx for exp an z inputs used in viz

    sig_gctx_val = out_directory+ "/z_scores.gct"
    # if args.zscore_gct :
    #     sig_gctx_val = args.zscore_gct


    # run eVIP_corr.py
    # print('calculating correlations...')
    run_corr = eVIP_corr.run_main(input=infile,zscore_gct=zscore_gct,
                                    out_dir= out_directory)

    # print('comparing...')
    run_compare = eVIP_compare.run_main(sig_info=sig_info,
                    gctx = out_directory+"/spearman_rank_matrix.gct",
                    allele_col = args.allele_col, o= out_directory+"/compare",
                    r = args.r, c = args.c, i = None,
                    conn_null = args.conn_null, ie_col = args.ie_col,
                    ie_filter = args.ie_filter, num_reps = args.num_reps,
                    cell_id = args.cell_id, plate_id = args.plate_id)

    # print('predicting...')
    run_predict = eVIP_predict.run_main(i= out_directory+"/compare.txt",
                o= out_directory+"/predict", conn_thresh=args.conn_thresh,
                mut_wt_rep_thresh=args.mut_wt_rep_thresh,
                mut_wt_rep_rank_diff=args.mut_wt_rep_rank_diff,
                disting_thresh=args.disting_thresh,
                use_c_pval=args.use_c_pval,
                cond_median_max_diff_thresh=args.cond_max_diff_thresh)


    if not args.sparkler_off:
        # print "making sparkler plots..."
        run_sparkler = eVIP_sparkler.eVIP_run_main(pred_file = out_directory+"/predict.txt",
                        ref_allele_mode=args.ref_allele_mode,y_thresh = args.y_thresh ,
                        x_thresh = args.x_thresh,use_c_pval= args.use_c_pval,
                        annotate=args.annotate, by_gene_color= args.by_gene_color,
                        pdf= args.pdf,xmin= args.xmin, xmax = args.xmax,
                        ymin = args.ymin, ymax = args.ymax,
                        out_dir = out_directory+"/sparkler_plots")

    if not args.viz_off:
        # print "making visualizations..."
        if args.conn_null:
            null_conn = args.conn_null
        else:
            null_conn = out_directory + "/compare_conn_null.txt"

        run_viz = eVIP_viz.eVIP_run_main(pred_file= out_directory+"/predict.txt",
                sig_info = args.sig_info,
                gctx=out_directory+"/spearman_rank_matrix.gct",
                sig_gctx = sig_gctx_val, ref_allele_mode = args.ref_allele_mode,
                null_conn = null_conn,out_dir = out_directory+"/viz",
                ymin = args.viz_ymin, ymax= args.viz_ymax,
                allele_col = args.allele_col, use_c_pval = args.use_c_pval,
                pdf = args.pdf, cell_id = args.cell_id, plate_id = args.plate_id,
                corr_val_str= args.corr_val)


#################
# END FUNCTIONS #
#################

if __name__ == "__main__": main()
