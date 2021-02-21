#!/usr/bin/python

#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# <Script name>
# Author: Angela Brooks
# Program Completion Date:
# Description:
# Modification Date(s):
# Copyright (c) 2011, Angela Brooks. anbrooks@gmail.com
# All rights reserved.

import sys
import optparse
import os
import pdb
import csv

###########
# CLASSES #
###########
class OptionParser(optparse.OptionParser):
    """
    Adding a method for required arguments.
    Taken from:
    http://www.python.org/doc/2.3/lib/optparse-extending-examples.html
    """
    def check_required(self, opt):
        option = self.get_option(opt)

        # Assumes the option's 'default' is set to None!
        if getattr(self.values, option.dest) is None:
            print "%s option not supplied" % option
            self.print_help()
            sys.exit(1)

###############
# END CLASSES #
###############

########
# MAIN #
########
def main():

    opt_parser = OptionParser()

    # Add Options. Required options should have default=None
    opt_parser.add_option("-i",
                          dest="input_table",
                          type="string",
                          help="""Input table with mutation impact comparison
                                  pvals""",
                          default=None)
    opt_parser.add_option("-o",
                          dest="output_table",
                          type="string",
                          help="""Output table with mutation impact predictions
                                  based on different methods""",
                          default=None)
    opt_parser.add_option("--conn_thresh",
                          dest="conn_thresh",
                          type="float",
                          help="P-value threshould for connectivity vs null.",
                          default=0.05)
    opt_parser.add_option("--mut_wt_rep_thresh",
                          dest="mut_wt_thresh",
                          type="float",
                          help="""P-value threshould for comparison of WT and mut
                                  robustness""",
                          default=0.05)
    opt_parser.add_option("--mut_wt_rep_rank_diff",
                          dest="mut_wt_rep_diff",
                          type="float",
                          help="""The minimum difference in median rankpoint
                                  between WT and mut to consider a difference.
                                  DEF=0""" ,
                          default=0)
    opt_parser.add_option("--disting_thresh",
                          dest="disting_thresh",
                          type="float",
                          help="""P-value threshould that tests if mut and wt reps
                                  are indistinguishable from each other""",
                          default=0.05)
    opt_parser.add_option("--use_c_pval",
                          dest="use_c_pval",
                          action="store_true",
                          help="Will use corrected p-value instead of raw p-val",
                          default=False)
    opt_parser.add_option("--cond_median_max_diff_thresh",
                          dest="cond_median_max_diff_thresh",
                          help="""Threshold for maximum difference between
                          condition correlation medians when determining if
                           variant is not neutral. Default = 0.2""",
                          default=0.2)

    (options, args) = opt_parser.parse_args()

    # validate the command line arguments
    opt_parser.check_required("-i")
    opt_parser.check_required("-o")
    opt_parser.check_required("--mut_wt_rep_thresh")
    opt_parser.check_required("--disting_thresh")


    run_main(i=options.input_table, o= options.output_table,
            conn_thresh=options.conn_thresh,
            mut_wt_rep_thresh=options.mut_wt_thresh,
            mut_wt_rep_rank_diff=options.mut_wt_rep_diff,
            disting_thresh=options.disting_thresh,
            use_c_pval=options.use_c_pval,
            cond_median_max_diff_thresh=options.cond_median_max_diff_thresh)


def run_main(i=None, o= None, conn_thresh=None, mut_wt_rep_thresh=None,
             mut_wt_rep_rank_diff=None, disting_thresh=None,
            use_c_pval=None, cond_median_max_diff_thresh=None):

    #setting default values
    c_thresh = float(conn_thresh) if conn_thresh != None else float(.1)
    mut_wt_thresh = float(mut_wt_rep_thresh) if mut_wt_rep_thresh != None else float(.1)
    mut_wt_rep_diff = float(mut_wt_rep_rank_diff) if mut_wt_rep_rank_diff != None else float(0)
    disting_thresh = float(disting_thresh) if disting_thresh != None else float(.1)
    cond_median_max_diff_thresh = float(cond_median_max_diff_thresh) if cond_median_max_diff_thresh != None else float(.2)

    input_table = open(i)
    output = open(o+".txt", "w")

    file_reader = csv.DictReader(input_table, delimiter="\t")

    column_headers = ["gene",
                      "mut",
                      "mut_rep",
                      "wt_rep",
                      "mut_wt_connectivity",
                      "wt",
                      "cell_line",
                      "mut_wt_rep_pval",
                      "mut_wt_conn_null_pval",
                      "wt_mut_rep_vs_wt_mut_conn_pval",
                      "kruskal_diff",
                      "mut_wt_rep_c_pval",
                      "mut_wt_conn_null_c_pval",
                      "wt_mut_rep_vs_wt_mut_conn_c_pval",
                      "prediction"]


    file_writer = csv.DictWriter(output, delimiter="\t",
                                 fieldnames=column_headers)
    file_writer.writeheader()


    for row in file_reader:
        if use_c_pval:
            prediction = get_prediction_6(float(row["wt_rep"]),
                                         float(row["mut_rep"]),
                                         float(row["mut_wt_rep_c_pval"]),
                                         float(row["mut_wt_connectivity"]),
                                         float(row["mut_wt_conn_null_c_pval"]),
                                         float(row["wt_mut_rep_vs_wt_mut_conn_c_pval"]),
                                         mut_wt_thresh,
                                         mut_wt_rep_diff,
                                         c_thresh,
                                         disting_thresh,
                                         cond_median_max_diff_thresh,
                                         float(row["kruskal_diff"]))
            row["prediction"]= prediction

        else:
            prediction = get_prediction_6(float(row["wt_rep"]),
                                         float(row["mut_rep"]),
                                         float(row["mut_wt_rep_pval"]),
                                         float(row["mut_wt_connectivity"]),
                                         float(row["mut_wt_conn_null_pval"]),
                                         float(row["wt_mut_rep_vs_wt_mut_conn_pval"]),
                                         mut_wt_thresh,
                                         mut_wt_rep_diff,
                                         c_thresh,
                                         disting_thresh,
                                         cond_median_max_diff_thresh,
                                         float(row["kruskal_diff"]))

            row["prediction"]= prediction

        file_writer.writerow(row)

    input_table.close()
    output.close()

    # sys.exit(0)


############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def get_prediction_6(wt_rep, mut_rep, mut_wt_rep_pval,
                     mut_wt_conn, mut_wt_conn_pval, disting_pval,
                     mut_wt_thresh, mut_wt_rep_diff, c_thresh, disting_thresh,
                     cond_median_max_diff_thresh,cond_median_max_diff):

    if disting_pval < disting_thresh:
        if max_diff(wt_rep, mut_rep, mut_wt_conn) < mut_wt_rep_diff:
            if mut_wt_conn_pval < c_thresh:
                return "Neutral"
            else:
                return "Error"

#       if mut_wt_conn_pval < c_thresh:
#           if mut_wt_conn < conn_null_med:
#               return "DOM-NEG"

        #is diff between the max and min of the 3 comparisons from kruskal> threshold
        if cond_median_max_diff > cond_median_max_diff_thresh:

            if mut_wt_rep_pval < mut_wt_thresh:
                if wt_rep < mut_rep:
                    if mut_rep - wt_rep >= mut_wt_rep_diff:
                        return "GOF"
                    else:
                        return "COF"
                elif wt_rep > mut_rep:
                    if wt_rep - mut_rep >= mut_wt_rep_diff:
                        return "LOF"
                    else:
                        return "COF"
                else:
                    return "COF"
            else:
                return "COF"

    if mut_wt_conn_pval < c_thresh:
        return "Neutral"

    return "NI"

def max_diff(wt_rep, mut_rep, mut_wt_conn):
    max_diff = abs(wt_rep - mut_rep)

    wt_conn_diff = abs(wt_rep - mut_wt_conn)
    if wt_conn_diff > max_diff:
        max_diff = wt_conn_diff

    mut_conn_diff = abs(mut_rep - mut_wt_conn)
    if mut_conn_diff > max_diff:
        max_diff = mut_conn_diff

    return max_diff
#################
# END FUNCTIONS #
#################
if __name__ == "__main__": main()
