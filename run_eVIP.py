#!/usr/bin/python
import sys
import argparse
import os
import errno
import csv
import itertools
import rpy2.robjects as robjects
import json

#importing eVIP
from bin import eVIP_corr
from bin import eVIP_predict
from bin import eVIP_sparkler
from bin import eVIP_viz
from bin import eVIP_compare

########
# MAIN #
########
def main(infile=None, zscore_gct = None, out_directory=None, sig_info =None,
            c=None, r=None, num_reps=None,ie_filter=None,ie_col=None, i=None,
            allele_col=None, conn_null=None, conn_thresh=None,
            mut_wt_rep_rank_diff=None, use_c_pval=None, cell_id=None,
            plate_id=None, ref_allele_mode=None,x_thresh=None, y_thresh=None,
            annotate=None, by_gene_color=None, pdf=None, xmin=None,xmax=None,
            ymin=None, ymax=None, viz_ymin=None, viz_ymax=None, corr_val=None):

    parser = argparse.ArgumentParser()

    parser.add_argument("-l1000", action ="store_true",
                        help = """Using l1000 data as input. Must use -zscore_gct""")


    #from corr
    parser.add_argument("--infile",
                        help="Input txt file (filtered and log transformed data).")
    parser.add_argument("-zscore_gct",
                        help="Zscore input gct file (use instead of --infile)")
    parser.add_argument("-out_directory",required=True,
                        help="Path to directory for eVIP output files")

    #from compare
    parser.add_argument("-sig_info",required=True,
                        help = """sig info file with gene information and distil
                        information""")
    parser.add_argument("-c",required=True,
                        help = """.grp file containing allele names of control
                        perturbations.If this file is given, a null will be
                        calculated from these""")
    parser.add_argument("-r", required=True,
                        help = """File explicitly indicating which comparisons
                        to do. Assumes the file has a header. The first column
                        is the reference allele and second column is test allele.""")
    parser.add_argument("-num_reps",required=True,
                        help = "Number of replicates expected for each allele. DEF=3")
    parser.add_argument("-ie_filter",
                        help = """L1000 threshold for infection efficiency. Any
                        wildtype or mutant alleles having an ie below this threshold,
                        will be removed""")
    parser.add_argument("-ie_col",
                        help = """Name of the column in the sig_info file with
                        infection efficiency information. DEF=x_ie_a549""")
    parser.add_argument("-i",
                        help = "Number of iterations to run. DEF=1000")
    parser.add_argument("-allele_col",
                        help = """Column name in sig_info file that indicates
                        the allele names. DEF=x_mutation_status""")
    parser.add_argument("-conn_null",
                        help = """Optional file containing connectivity null
                        values from a previous run. Should end in _conn_null.txt""")

    #from predict
    parser.add_argument("-conn_thresh",
                        help = """P-value threshold for connectivity vs null.
                                DEFAULT=0.05""")
    parser.add_argument("-mut_wt_rep_thresh",
                        help = """P-value threshold for comparison of WT and
                        mut robustness. DEFAULT=0.05""")
    parser.add_argument("-disting_thresh",
                        help = """P-value threshold that tests if mut and wt reps
                        are indistinguishable from each other.DEFAULT=0.05""")
    parser.add_argument("-mut_wt_rep_rank_diff",
                        help = """The minimum difference in median rankpoint WT
                        and mut to consider a difference. DEF=0""")
    parser.add_argument("-use_c_pval", action ="store_true",
                        help = """Will use corrected p-value instead of raw p-val""")
    parser.add_argument("-cell_id", help = """Optional: Will only look at
                        signatures from this cell line. Helps to filter
                        sig_info file.""")
    parser.add_argument("-plate_id", help = """Optional: Will only look at
                        signatures from this plate""")

    #from sparkler
    parser.add_argument("-ref_allele_mode", action ="store_true",
                        help = """Sparkler+Viz: Instead of organizing plots by
                        gene, will use the wt column to determine what are the
                        reference alleles.""" )
    parser.add_argument("-x_thresh" ,
                        help = "Sparkler: Threshold of significance")
    parser.add_argument("-y_thresh",
                        help = "Sparkler: Threshold of impact direction")
    parser.add_argument("-annotate", action ="store_true",
                        help = "Sparkler: Will add allele labels to points.")
    parser.add_argument("-by_gene_color",
                        help = """Sparkler: File containing labels and
                        colors for gene-centric plot.""")
    parser.add_argument("-pdf",
                        help = """Sparkler + Viz: Will print plots in pdf
                        format instead of png.""")
    parser.add_argument("-xmin",
                        help = """Sparkler: Min value of x-axis. DEF=0""")
    parser.add_argument("-xmax",
                        help = """Sparkler: Max value of x-axis. DEF=4""")
    parser.add_argument("-ymin",
                        help = "Sparkler: Min value of y-axis. DEF=-3")
    parser.add_argument("-ymax",
                        help = "Sparkler: Min value of y-axis. DEF=3")
    #from viz
    parser.add_argument("-viz_ymin",
                        help = "Viz: Minimum y-value of rep value. DEF=-100")
    parser.add_argument("-viz_ymax",
                        help = "Viz: Maximum y-value of rep value. DEF=100")
    parser.add_argument("-corr_val",
                        help = """Viz: String used to label the correlation
                        value. DEF= 'row median rankpoints' """)

    parser.add_argument("-sparkler_off", action ="store_true",
                        help = "Will not perform eVIP sparkler step")
    parser.add_argument("-viz_off", action ="store_true",
                        help = "Will not perform eVIP viz step")


    global args
    args = parser.parse_args()

    #make eVIP output directory
    global out_dir
    out_dir = args.out_directory
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)


        eVIP_dir = out_dir + "/eVIP_out"
        if not os.path.exists(eVIP_dir):
            os.makedirs(eVIP_dir)





    if args.l1000:

        run_eVIP(args.infile, args.zscore_gct, eVIP_dir, args.sig_info, args.c,
                args.r, args.num_reps,args.ie_filter, args.ie_col, args.i,
                args.allele_col, args.conn_null, args.conn_thresh,
                args.mut_wt_rep_rank_diff, args.use_c_pval, args.cell_id, args.plate_id, args.ref_allele_mode,
                 args.x_thresh, args.y_thresh, args.annotate, args.by_gene_color, args.pdf, args.xmin,
                 args.xmax, args.ymin, args.ymax, args.viz_ymin, args.viz_ymax, args.corr_val)

############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def run_eVIP(infile=None, zscore_gct = None, out_directory=None,
            sig_info =None, c=None, r=None, num_reps=None,
         ie_filter=None,ie_col=None, i=None, allele_col=None, conn_null=None, conn_thresh=None,
         mut_wt_rep_rank_diff=None, use_c_pval=None, cell_id=None, plate_id=None, ref_allele_mode=None,
         x_thresh=None, y_thresh=None, annotate=None, by_gene_color=None, pdf=None, xmin=None,
         xmax=None, ymin=None, ymax=None, viz_ymin=None, viz_ymax=None, corr_val=None):



    #different sig_gctx for exp an z inputs used in viz
    if args.infile :
        sig_gctx_val = out_directory+ "/z_scores.gct"
    if args.zscore_gct :
        sig_gctx_val = args.zscore_gct


    # run eVIP_corr.py
    print('calculating correlations...')
    run_corr = eVIP_corr.run_main(input=infile,zscore_gct=zscore_gct, out_dir= out_directory)

    print('comparing...')
    run_compare = eVIP_compare.run_main(sig_info=sig_info, gctx = out_directory+"/spearman_rank_matrix.gct",
                allele_col = args.allele_col, o= out_directory+"/compare", r = args.r,
             c = args.c, i = args.i, conn_null = args.conn_null, ie_col = args.ie_col,
             ie_filter = args.ie_filter, num_reps = args.num_reps, cell_id = args.cell_id, plate_id = args.plate_id)

    print('predicting...')
    run_predict = eVIP_predict.run_main(i= out_directory+"/compare.txt", o= out_directory+"/predict", conn_thresh=args.conn_thresh,
                mut_wt_rep_thresh=args.mut_wt_rep_thresh, mut_wt_rep_rank_diff=args.mut_wt_rep_rank_diff,
                disting_thresh=args.disting_thresh, use_c_pval=args.use_c_pval)


    if not args.sparkler_off:
        print "making sparkler plots..."
        run_sparkler = eVIP_sparkler.eVIP_run_main(pred_file = out_directory+"/predict.txt", ref_allele_mode=args.ref_allele_mode,
                y_thresh = args.y_thresh , x_thresh = args.x_thresh,
                use_c_pval= args.use_c_pval,annotate=args.annotate, by_gene_color= args.by_gene_color, pdf= args.pdf,
                xmin= args.xmin, xmax = args.xmax, ymin = args.ymin, ymax = args.ymax, out_dir = out_directory+"/sparkler_plots")

    if not args.viz_off:
        print "making visualizations..."
        if args.conn_null:
            null_conn = args.conn_null
        else:
            null_conn = out_directory + "/compare_conn_null.txt"

        run_viz = eVIP_viz.eVIP_run_main(pred_file= out_directory+"/predict.txt", sig_info = args.sig_info, gctx=out_directory+"/spearman_rank_matrix.gct",
                sig_gctx = sig_gctx_val, ref_allele_mode = args.ref_allele_mode, null_conn = null_conn,
                out_dir = out_directory+"/viz",ymin = args.viz_ymin, ymax= args.viz_ymax, allele_col = args.allele_col, use_c_pval = args.use_c_pval,
                 pdf = args.pdf, cell_id = args.cell_id, plate_id = args.plate_id, corr_val_str= args.corr_val)


#################
# END FUNCTIONS #
#################

if __name__ == "__main__":
    main()
