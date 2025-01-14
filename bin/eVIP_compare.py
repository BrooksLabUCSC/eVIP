#!/usr/bin/python


#!/broad/software/free/Linux/redhat_5_x86_64/pkgs/python_2.5.4/bin/python
# mutation_impact
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
import math
import random
import numpy
import rpy2.robjects as robjects
import cmapPy.pandasGEXpress.parse as gct
import cmapPy.pandasGEXpress.subset as grp
from cmapPy.pandasGEXpress.parse import parse

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

#############
# CONSTANTS #
#############
NUM_ITERATIONS = 1000
NUM_REPS = 3

#LOG10_ZERO = 10.0
LOG10_ZERO = 35.0

# DEF_IE_COL = "x_ie_a549"
# DEF_ALLELE_COL = "x_mutation_status"

#################
# END CONSTANTS #
#################


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
    opt_parser.add_option("--sig_info",
                          dest="sig_info",
                          type="string",
                          help="""sig info file with gene information and distil
                                  information""",
                          default=None)
    opt_parser.add_option("--gctx",
                          dest="gctx",
                          type="string",
                          help="""GCTX file of pairwise similarity. Designed to
                                  use connectivity score, but any measurement of
                                  similarity can be used (e.g., pearson
                                  correlation)""",
                          default=None)
    opt_parser.add_option("--allele_col",
                          dest="allele_col",
                          type="string",
                          help="""Column name that indicates the allele names.
                                  DEF= "allele" """ ,
                          default="allele")
    opt_parser.add_option("-o",
                          dest="output_file_prefix",
                          type="string",
                          help="""Prefix of output files of mutation impact data.
                                  Includes figures of p-value distribution.""",
                          default=None)
    opt_parser.add_option("-r",
                          dest="reference_test_file",
                          type="string",
                          help="""File explicitly indicating which comparisons
                                  to do. Assumes the file has a header and it is
                                  ignored. The first column is the reference
                                  allele and second column is test allele. If
                                  this file is not given, then the reference
                                  alleles are assumed to be WT and inferred from
                                  the allele names.""",
                          default=None)
    opt_parser.add_option("-c",
                          dest="controls_file",
                          type="string",
                          help=""".grp file containing allele names of control
                                  perturbations. If this file is given, a null
                                  will be calculated from these""",
                          default=None)
    opt_parser.add_option("-i",
                          dest="num_iterations",
                          type="int",
                          help="Number of iterations to run. DEF=%d" % NUM_ITERATIONS,
                          default=NUM_ITERATIONS)
    opt_parser.add_option("--conn_null",
                          dest="conn_null_input",
                          type="string",
                          help="""Optional file containing connectvity null
                          values from a previous run. Should end
                          in _conn_null.txt""",
                          default=None)
    opt_parser.add_option("--ie_col",
                          dest="ie_col",
                          type="string",
                          help="""Name of the column with infection efficiency
                                  information.""" )
    opt_parser.add_option("--ie_filter",
                          dest="ie_filter",
                          type="float",
                          help="""Threshold for infection efficiency. Any wildtype
                                  or mutant alleles having an ie below this
                                  threshold, will be removed""",
                          default=None)
    opt_parser.add_option("--num_reps",
                          dest="num_reps",
                          type="int",
                          help="""Number of replicates expected for each allele.
                                  DEF=%d""" % NUM_REPS,
                          default=NUM_REPS)
    opt_parser.add_option("--cell_id",
                          dest="cell_id",
                          type="string",
                          help="""Optional: Will only look at signatures from
                           this cell line. Helps to filter sig_info file.""",
                          default=None)
    opt_parser.add_option("--plate_id",
                          dest="plate_id",
                          type="string",
                          help="""Optional: Will only look at signatures from
                                  this plate.""",
                          default=None)

    (options, args) = opt_parser.parse_args()

    # validate the command line arguments
    opt_parser.check_required("--sig_info")
    opt_parser.check_required("--gctx")
    opt_parser.check_required("-c")
    opt_parser.check_required("-o")



    run_main(sig_info=options.sig_info, gctx = options.gctx,
                allele_col = options.allele_col, o = options.output_file_prefix,
                r = options.reference_test_file, c = options.controls_file,
                i = options.num_iterations, conn_null = options.conn_null_input,
                ie_col = options.ie_col,  ie_filter = options.ie_filter,
                num_reps = options.num_reps, cell_id = options.cell_id,
                plate_id = options.plate_id)


def run_main(sig_info=None, gctx = None, allele_col = None, o = None, r = None,
             c = None, i = None, conn_null = None, ie_col = None,
             ie_filter = None, num_reps = None, cell_id = None, plate_id = None):


    #default values
    i = int(i) if i != None else int(1000)
    # ie_col = str(ie_col) if ie_col != None else str(x_ie_a549)
    ie_filter = float(ie_filter) if ie_filter != None else float(0.0)
    if ie_filter == float(0.0):
        ie_col = None
    num_reps = int(num_reps) if num_reps != None else int(3)


    sig_info_file = open(sig_info)
    output_file_prefix = open(o + ".txt", "w")

    # Output distribution files
    controls = grp.grp.read(c)

    reference_test_filename = r
    ref2test_allele = None
    if reference_test_filename:
        ref2test_allele = parseRefTestFile(reference_test_filename)

    if ref2test_allele == None:
        print("Error reading in comparisons file")
        sys.exit()

    this_gctx = parse(gctx)
    # this_gctx.read()

    num_iterations = int(i)
    num_reps = int(num_reps)

    conn_null_input = conn_null

    if conn_null_input:
        conn_nulls_from_input_str = grp.grp.read(conn_null_input)
        conn_nulls_from_input = map(float, conn_nulls_from_input_str)

    (allele2distil_id,
     allele2WT,
     allele2gene,
     allele2cell_id,
     WT_alleles) = parse_sig_info(sig_info_file,
                                  ref2test_allele,
                                  allele_col,
                                  ie_col, ie_filter,
                                  cell_id,
                                  plate_id)

    clean_controls = []
    for this_control in controls:
        if this_control in allele2distil_id:
            clean_controls.append(this_control)


    #calculates if no inputs
    replicate_null_dist, connectivity_null_dist = getNullDist(this_gctx,
                                                            allele2distil_id,
                                                            clean_controls,
                                                            num_iterations,
                                                            num_reps)



    #overwrites conn_null_dist if its an input
    if conn_null_input:
        connectivity_null_dist = conn_nulls_from_input

    if not conn_null:
        conn_null_dist_out = open(o + "_conn_null.txt", "w")
        for x in connectivity_null_dist:
            conn_null_dist_out.write("%f\n" % x)
        conn_null_dist_out.close()



    WT_dict, wt_rep_pvals, wt_ordered = buildWT_dict(this_gctx,
                                                    allele2distil_id, WT_alleles,
                                                    replicate_null_dist, num_reps)

    # Print header to output file
    output_file_prefix.write("gene\tmut\tmut_rep\twt_rep\tmut_wt_connectivity\t")
    output_file_prefix.write("wt\tcell_line\t")
    output_file_prefix.write("mut_wt_rep_pval\tmut_wt_conn_null_pval\twt_mut_rep_vs_wt_mut_conn_pval\tkruskal_diff\t")
    output_file_prefix.write("mut_wt_rep_c_pval\tmut_wt_conn_null_c_pval\twt_mut_rep_vs_wt_mut_conn_c_pval\n")

    mut_rep_pvals = []
    mut_wt_rep_pvals = []
    mut_wt_conn_pvals = []
    mut_wt_rep_vs_wt_mut_conn_pvals = []

    outlines = []

    # Build comparison
    for allele in allele2WT:

        # Don't calculate for the WT allele
        if allele == allele2WT[allele]:
            continue

        mut_rankpt, mut_rankpt_dist = getSelfConnectivity(this_gctx,
                                                    allele2distil_id[allele],
                                                    num_reps)

        self_pval = getPairwiseComparisons(mut_rankpt_dist,
                                           replicate_null_dist)
        mut_rep_pvals.append(self_pval)

        mut_wt_conn_rankpt, mut_wt_conn_dist = getConnectivity(this_gctx,
                                                    allele2distil_id[allele],
                                                    allele2distil_id[allele2WT[allele]],
                                                    num_reps)

        conn_pval = getPairwiseComparisons(mut_wt_conn_dist,
                                           connectivity_null_dist)
        mut_wt_conn_pvals.append(conn_pval)

        mut_wt_rep_pval = getPairwiseComparisons(mut_rankpt_dist,
                                                 WT_dict[allele2WT[allele]]["wt_rep_dist"])
        mut_wt_rep_pvals.append(mut_wt_rep_pval)


        wt_mut_rep_vs_wt_mut_conn_pval = getKruskal(WT_dict[allele2WT[allele]]["wt_rep_dist"],
                                                    mut_rankpt_dist,
                                                    mut_wt_conn_dist)
        mut_wt_rep_vs_wt_mut_conn_pvals.append(wt_mut_rep_vs_wt_mut_conn_pval)


        medians = []
        medians.append(median(WT_dict[allele2WT[allele]]["wt_rep_dist"]))
        medians.append(median(mut_rankpt_dist))
        medians.append(median(mut_wt_conn_dist))

        median_diff = max(medians)-min(medians)

        out_elems = [allele2gene[allele],
                     allele,
                     "%f" % mut_rankpt,
                     "%f" % WT_dict[allele2WT[allele]]["wt_rep"],
                     "%f" % mut_wt_conn_rankpt,
                     allele2WT[allele],
                     allele2cell_id[allele],
                     "%f" % mut_wt_rep_pval,
                     "%f" % conn_pval,
                     "%f" % wt_mut_rep_vs_wt_mut_conn_pval,
                     "%f" % median_diff]
        outline = "\t".join(out_elems)
        outlines.append(outline)

    # Calculate corrected pvalues
    mut_wt_rep_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(mut_wt_rep_pvals), "BH")
    mut_wt_conn_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(mut_wt_conn_pvals), "BH")
    mut_wt_rep_vs_wt_mut_conn_c_pvals = robjects.r['p.adjust'](robjects.FloatVector(mut_wt_rep_vs_wt_mut_conn_pvals),
                                                               "BH")

    # Write to file
    num_lines = len(outlines)
    for i in range(num_lines):
        this_outline = outlines[i]

        this_outline += "\t%f\t" % mut_wt_rep_c_pvals[i]
        this_outline += "%f\t" % mut_wt_conn_c_pvals[i]
        this_outline += "%f\n" % mut_wt_rep_vs_wt_mut_conn_c_pvals[i]


        output_file_prefix.write(this_outline)


############
# END_MAIN #
############

#############
# FUNCTIONS #
#############
def median(lst):
    n = len(lst)
    s = sorted(lst)
    return (sum(s[n//2-1:n//2+1])/2.0, s[n//2])[n % 2] if n else None

def buildWT_dict(this_gctx, allele2distil_id, WT_alleles,
                    replicate_null_dist, num_reps):
    """
    {WT_allele:{"wt_rep": med_wt_rep,
                "wt_rep_dist":[]
                "wt_rep_pval": p_val vs null}
    """
    WT_dict = {}
    wt_rep_pvals = []
    wt_allele_ordered = []
    for allele in WT_alleles:
        WT_dict[allele] = {}
        wt_rep_rankpt, rep_rankpts = getSelfConnectivity(this_gctx,
                                                        allele2distil_id[allele],
                                                        num_reps)

        WT_dict[allele]["wt_rep"] = wt_rep_rankpt
        WT_dict[allele]["wt_rep_dist"] = rep_rankpts
        wt_rep_pval = getPairwiseComparisons(rep_rankpts,
                                             replicate_null_dist)
        WT_dict[allele]["wt_rep_pval"] = wt_rep_pval
        wt_rep_pvals.append(wt_rep_pval)
        wt_allele_ordered.append(allele)

    return WT_dict, wt_rep_pvals, wt_allele_ordered

def formatDir(i_dir):
    i_dir = os.path.realpath(i_dir)
    if i_dir.endswith("/"):
        i_dir = i_dir.rstrip("/")
    return i_dir

def formatLine(line):
    line = line.replace("\r","")
    line = line.replace("\n","")
    return line

def getKruskal(wt_rankpt_dist, mut_rankpt_dist, mut_wt_conn_dist):
    return robjects.r["kruskal.test"](robjects.ListVector(
                                    {'a':robjects.FloatVector(wt_rankpt_dist),
                                    'b':robjects.FloatVector(mut_rankpt_dist),
                                    'c':robjects.FloatVector(mut_wt_conn_dist)}))[2][0]


def getSelfConnectivity(this_gctx, distil_ids, num_reps):
    """
    returns
    median rankpoint
    Distribution of values
    """
    row_medians = []
    for i in range(num_reps):
        rank_pts = []
        for j in range(num_reps):
            if i == j:
                continue
            rank_pts.append(float(this_gctx.data_df[distil_ids[i]]
                                                 [distil_ids[j]]))

        row_medians.append(numpy.percentile(rank_pts, 50))

    return numpy.percentile(row_medians, 50), row_medians


#   rep_rankpts = []
#   for i in range(num_reps):
#       for j in range(i+1, num_reps):
#           rep_rankpts.append(float(this_gctx.data_df[distil_ids[i]]
#                                                   [distil_ids[j]]))
#
#   return numpy.percentile(rep_rankpts, 50), rep_rankpts

def getConnectivity(this_gctx, distil_ids1, distil_ids2, num_reps):
    """
    Returns
    median connectivity
    Distribution of values
    """
    row_column_medians = []

    col_rankpts = [[] for n in range(num_reps)]
    for i in range(num_reps):
        row_rankpts = []
        for j in range(num_reps):
            val = float(this_gctx.data_df[distil_ids1[i]][distil_ids2[j]])
            row_rankpts.append(val)
            col_rankpts[j].append(val)

        row_column_medians.append(numpy.percentile(row_rankpts, 50))

    for n in range(num_reps):
        row_column_medians.append(numpy.percentile(col_rankpts[n], 50))

    return numpy.percentile(row_column_medians, 50), row_column_medians

#   conn_rankpts = []
#   for i in range(num_reps):
#       for j in range(num_reps):
#           conn_rankpts.append(float(this_gctx.data_df[distil_ids1[i]]
#                                                    [distil_ids2[j]]))

#   return numpy.percentile(conn_rankpts, 50), conn_rankpts

def getLog(vals):
    out_vals = []
    for val in vals:
        if val == 0:
            out_vals.append(LOG10_ZERO)
            continue
        out_vals.append(-math.log10(val))

    return out_vals

def getPairwiseComparisons(dist1, dist2):

    return robjects.r["wilcox.test"](robjects.FloatVector(dist1),
                                     robjects.FloatVector(dist2))[2][0]

def getNullDist(this_gctx, allele2distil_id, controls,
                num_iterations, num_reps):
    """
    Returns
    replicate_null_dist,
    connectivity_null_dist
    """
    rep_null_dist = []
    connectivity_null_dist = []

    allele_list = allele2distil_id.keys()
    #get combination of each

    all_combos = [(x,y) for x in controls for y in allele_list]

    for c in all_combos:

        random_control = c[0]
        random_allele = c[1]

        if random_control == random_allele:
            continue

        control_distil_ids_set = set([random.choice(allele2distil_id[random_control])])
        while len(control_distil_ids_set) < num_reps:
            this_control = random.choice(controls)
            this_distil = random.choice(allele2distil_id[this_control])
            if this_distil in control_distil_ids_set:
                continue
            control_distil_ids_set.add(this_distil)

        control_distil_ids = list(control_distil_ids_set)

        # Introspect similarity is median of row similarities
        rank_pts = []
        for i in range(1, num_reps):

            rank_pts.append(float(this_gctx.data_df[control_distil_ids[0]][control_distil_ids[i]]))

        rep_null_dist.append(numpy.percentile(rank_pts, 50))

        # Get connectivity null
        rank_pts = []
        for i in range(num_reps):
            rank_pts.append(float(this_gctx.data_df[allele2distil_id[random_control][0]]
                                                 [allele2distil_id[random_allele][i]]))

        connectivity_null_dist.append(numpy.percentile(rank_pts, 50))


    return rep_null_dist, connectivity_null_dist


def hasLowIE(ie_string, ie_thresh):
    ie_elems = map(float, ie_string.split("|"))
    for ie in ie_elems:
        if ie < ie_thresh:
            return True
    return False


def parseRefTestFile(reference_test_filename):

    ref2test = {}

    ref_test_file = open(reference_test_filename)
    ctr = 1
    for line in ref_test_file:
        if ctr == 1:
            ctr += 1
            continue

        line = formatLine(line)
        lineList = line.split()

        updateDictOfLists(ref2test, lineList[0], lineList[1])

    return ref2test

def parse_sig_info(sig_info_file, ref2test_allele, allele_col, ie_col,
                    ie_filter = None, cell_id = None, plate_id = None):
    """
    Returns:
    allele2distil_id = {}
    allele2WT = {}
    allele2gene,
    allele2cell_id
    WT_alleles = []
    """
    allele2distil_id = {}
    gene2alleles = {}
    allele2WT = {}

    allele2gene = {}
    allele2cell_id = {}

    ie_col_idxs = []

    passed_alleles = set([])

    for line in sig_info_file:
        line = formatLine(line)
        lineList = line.split("\t")

        if "distil_id" in line:
            # sig_id_idx = lineList.index("sig_id")
            distil_idx = lineList.index("distil_id")
            gene_idx = lineList.index("pert_mfc_desc")
            allele_idx = lineList.index(allele_col)
            cell_idx = lineList.index("cell_id")

            if ie_filter:
                ie_cols = ie_col.split(",")
                for ie_col in ie_cols:
                    ie_col_idxs.append(lineList.index(ie_col))
            continue

        if cell_id:
            if cell_id != lineList[cell_idx]:
                continue

        # if plate_id:
        #     if plate_id not in lineList[sig_id_idx]:
        #         continue

        if ie_filter:
            for ie_idx in ie_col_idxs:
                low_ie_flag = hasLowIE(lineList[ie_idx], ie_filter)
                if low_ie_flag:
                    break
            if low_ie_flag:
                continue

        x_mutation_status = lineList[allele_idx]
        gene = lineList[gene_idx]

        if x_mutation_status in allele2distil_id:
            print "Warning: Duplicate allele present %s" % x_mutation_status
            continue

        distil_ids = lineList[distil_idx].split("|")

        for distil_id in distil_ids:
            updateDictOfLists(allele2distil_id, x_mutation_status, distil_id)

        updateDictOfSets(gene2alleles, gene, x_mutation_status)

        allele2gene[x_mutation_status] = gene
        allele2cell_id[x_mutation_status] = lineList[cell_idx]
        passed_alleles.add(x_mutation_status)


    if ref2test_allele:
        WT_alleles = set(ref2test_allele.keys())
        WT_alleles = WT_alleles & passed_alleles

        for ref in ref2test_allele:
            if ref not in passed_alleles:
                continue
            for test in ref2test_allele[ref]:
                if test not in passed_alleles:
                    continue
                allele2WT[test] = ref
    else: # Need to infer WT reference
        # Now find WT orf
        WT_alleles = set([])

        for gene in gene2alleles:
            # Find WT allele
            this_WT = ""
            if gene + "_WT.c" in gene2alleles[gene]:
                this_WT = gene + "_WT.c"
            elif gene + "_WT.c.2" in gene2alleles[gene]:
                this_WT = gene + "_WT.c.2"
            elif gene + "_WT.o" in gene2alleles[gene]:
                this_WT = gene + "_WT.o"
            elif gene + "_WT" in gene2alleles[gene]:
                this_WT = gene + "_WT"
            elif gene + "_WT.1" in gene2alleles[gene]:
                this_WT = gene + "_WT.1"
            elif gene + "_WT.2" in gene2alleles[gene]:
                this_WT = gene + "_WT.2"
            else:
                print "Warning: Can't find WT orf for %s" % gene
                continue

            WT_alleles.add(this_WT)

            # Match up alleles to WT
            for allele in gene2alleles[gene]:
                allele2WT[allele] = this_WT

    return allele2distil_id, allele2WT, allele2gene, allele2cell_id, list(WT_alleles)

def updateDictOfLists(d, key, item):
    """
    """
    try:
        d[key].append(item)
    except KeyError:
        d[key] = [item]

def updateDictOfSets(d, key, item):
    """
    Similar functionality as updateDictOfLists
    """
    try:
        d[key].add(item)
    except KeyError:
        d[key] = set([item])

#################
# END FUNCTIONS #
#################
if __name__ == "__main__": main()
