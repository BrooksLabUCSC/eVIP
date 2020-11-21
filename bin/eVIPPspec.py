########################################################################
# author: Alexis M. Thornton
#
########################################################################

import json
import os
import errno
import csv
import itertools
import rpy2.robjects as robjects

from bin import eVIP_corr
from bin import eVIP_predict
from bin import eVIP_sparkler
from bin import eVIP_viz
from bin import eVIP_compare
from bin import eVIPP_sparkler

def main(out_dir,JSON_file,min_genes,infile,sig_info, c, r, num_reps,
ie_filter, ie_col, i, allele_col, conn_null, conn_thresh,
mut_wt_rep_rank_diff, use_c_pval, cell_id, plate_id, ref_allele_mode,
x_thresh, y_thresh, annotate, by_gene_color, pdf, xmin,
xmax, ymin, ymax, viz_ymin, viz_ymax, corr_val,mut_wt_rep_thresh,disting_thresh,sparkler_off,viz_off,cond_median_max_diff_thresh):


    #getting used pathways (removing pathways in JSON that there isn't data for)
    pway_dict, used_pathways = JSON_pway(JSON_file,out_dir,infile,min_genes)

    if len(used_pathways) > 0:

        #running eVIP for each pathway
        JSON_eVIP(used_pathways,out_dir,infile,sig_info, c, r, num_reps,
        ie_filter, ie_col, i, allele_col, conn_null, conn_thresh,
        mut_wt_rep_rank_diff, use_c_pval, cell_id, plate_id, ref_allele_mode,
        x_thresh, y_thresh, annotate, by_gene_color, pdf, xmin,
        xmax, ymin, ymax, viz_ymin, viz_ymax, corr_val,mut_wt_rep_thresh,disting_thresh,sparkler_off,viz_off,cond_median_max_diff_thresh)

        summary_file = open(out_dir + "/eVIPP_summary.txt", "w")
        summarize(used_pathways, summary_file,out_dir)

        summary_eVIPP_vals = open(out_dir + "/eVIPP_combined_predict_files.txt", "w")
        summarize_predict_files(used_pathways, summary_eVIPP_vals,out_dir)

        #eVIPP sparkler
        # print "\n0"

        # print "Making allele pathway sparkler plots..."
        run_eVIPP_sparkler = eVIPP_sparkler.main(out_dir + "/eVIPP_combined_predict_files.txt", ref_allele_mode, y_thresh, x_thresh, use_c_pval,
                 annotate, by_gene_color, pdf, xmin, xmax, ymin, ymax,
                 out_dir + "/eVIPP_sparkler_plots")

    else:
        print("No pathways to test...")

if __name__ == "__main__":
    main()

def JSON_eVIP(used_pathways,out_dir,infile,sig_info, c, r, num_reps,
ie_filter, ie_col, i, allele_col, conn_null, conn_thresh,
mut_wt_rep_rank_diff, use_c_pval, cell_id, plate_id, ref_allele_mode,
x_thresh, y_thresh, annotate, by_gene_color, pdf, xmin,
xmax, ymin, ymax, viz_ymin, viz_ymax, corr_val,mut_wt_rep_thresh,disting_thresh,sparkler_off,viz_off,cond_median_max_diff_thresh):

    pways_mut_wt_rep_pvals_from_compare = []
    pways_mut_wt_conn_null_pvals_from_compare = []
    pways_wt_mut_rep_vs_wt_mut_conn_pvals_from_compare = []

    for pathway in used_pathways:
        pathway_out_file_txt = out_dir + "/" + pathway + "_eVIPP_outputs/" + pathway + "_expression.txt"
        z_pathway_out_file_txt = out_dir + "/" + pathway + "_eVIPP_outputs/" + pathway + "_zscores.gct"
        p_out_directory = out_dir + "/" + pathway + "_eVIPP_outputs"
        mkdir_p(p_out_directory)

        with open(out_dir + "/" + pathway + "_eVIPP_outputs/" + pathway + "_expression.txt","r") as pathway_out_file:

            if use_c_pval:
                # print "Running with multiple testing..."

                # running eVIP corr and compare to get pvals
                mut_wt_rep_pvals_from_compare, mut_wt_conn_null_pvals_from_compare, wt_mut_rep_vs_wt_mut_conn_pvals_from_compare,num_test_alleles \
                    = run_eVIP_multiple_testing_pt1(pathway_out_file_txt, None, p_out_directory, sig_info, c, r, num_reps,
                    ie_filter, ie_col, i, allele_col, conn_null, conn_thresh,
                    mut_wt_rep_rank_diff, use_c_pval, cell_id, plate_id, ref_allele_mode,
                    x_thresh, y_thresh, annotate, by_gene_color, pdf, xmin,
                    xmax, ymin, ymax, viz_ymin, viz_ymax, corr_val)


                # appending values from each pathway to one list
                pways_mut_wt_rep_pvals_from_compare.append(mut_wt_rep_pvals_from_compare)
                pways_mut_wt_conn_null_pvals_from_compare.append(mut_wt_conn_null_pvals_from_compare)
                pways_wt_mut_rep_vs_wt_mut_conn_pvals_from_compare.append(wt_mut_rep_vs_wt_mut_conn_pvals_from_compare)

                # if the pathway lists are done being collected
                if len(pways_mut_wt_rep_pvals_from_compare) == len(used_pathways):
                    eVIPP_mut_wt_rep_pvals = padj(pways_mut_wt_rep_pvals_from_compare, num_test_alleles, used_pathways)

                if len(pways_mut_wt_conn_null_pvals_from_compare) == len(used_pathways):
                    eVIPP_pways_mut_wt_conn_null_pvals = padj(pways_mut_wt_conn_null_pvals_from_compare,num_test_alleles, used_pathways)

                if len(pways_wt_mut_rep_vs_wt_mut_conn_pvals_from_compare) == len(used_pathways):
                    eVIPP_pways_wt_mut_rep_vs_wt_mut_conn_pvals = padj(pways_wt_mut_rep_vs_wt_mut_conn_pvals_from_compare,num_test_alleles, used_pathways)

                    #adding the new calculated eVIPP pathways to a new compare file
                    make_adj_compare_file(used_pathways, out_dir, eVIPP_mut_wt_rep_pvals,eVIPP_pways_mut_wt_conn_null_pvals,eVIPP_pways_wt_mut_rep_vs_wt_mut_conn_pvals)


                    for pathway in used_pathways:
                        # print "pathway running pt2 of eVIP"
                        # print pathway

                        pathway_out_file_txt = out_dir + "/" + pathway + "_expression.txt"
                        p_out_directory = out_dir + "/" + pathway + "_eVIPP_outputs"

                        #getting predictions and finisihin eVIP with the new corrected values
                        run_eVIP_multiple_testing_pt2(pathway_out_file_txt, None, p_out_directory, sig_info, c,
                                                      r, num_reps, ie_filter, ie_col, i, allele_col,
                                                      conn_null, conn_thresh,mut_wt_rep_rank_diff, use_c_pval,
                                                      cell_id,plate_id, ref_allele_mode,x_thresh, y_thresh,
                                                      annotate, by_gene_color,pdf, xmin,xmax, ymin,
                                                      ymax, viz_ymin, viz_ymax, corr_val,mut_wt_rep_thresh,disting_thresh,sparkler_off,viz_off,cond_median_max_diff_thresh)

                else:
                    run_eVIP(pathway_out_file_txt, None, p_out_directory, sig_info, c, r,
                             num_reps,
                             ie_filter, ie_col, i, allele_col, conn_null, conn_thresh,
                             mut_wt_rep_rank_diff, use_c_pval, cell_id, plate_id,
                             ref_allele_mode,
                             x_thresh, y_thresh, annotate, by_gene_color, pdf, xmin,
                             xmax, ymin, ymax, viz_ymin, viz_ymax, corr_val,mut_wt_rep_thresh,disting_thresh,sparkler_off,viz_off,cond_median_max_diff_thresh)


    return(used_pathways)

def JSON_pway(JSON_file,out_dir,infile,min_genes):
    pway_dict = {}
    with open(JSON_file, "rb") as JSON:
        old_dict = json.load(JSON)
        #creating new dict to remove pathways with None values (no genes were in the pathway)
        pway_dict = {k: v for k, v in old_dict.iteritems() if v is not None}

    used_pathways = []

    #for each pathway
    for pathway, genes in pway_dict.iteritems():
        print pathway

        pathway_out_file_txt = out_dir + "/" + pathway + "_eVIPP_outputs/" + pathway + "_expression.txt"
        p_out_directory = out_dir + "/" + pathway + "_eVIPP_outputs"
        mkdir_p(p_out_directory)

        pathway_genes = []
        for gene in genes:
            gene = str(gene)
            pathway_genes.append(gene)

        #running eVIP when the input is expression data
        with open(out_dir + "/" + pathway + "_eVIPP_outputs/" + pathway + "_expression.txt", "w+") as pathway_out_file, open(infile, "r") as all_data:
            matched_ids = []
            for line in all_data:
                if line.startswith("#gene_id"):
                    pathway_out_file.write(line)
                sline = line.split()

                for id in pathway_genes:
                    if sline[0]==id:
                        matched_ids.append(id)
                        pathway_out_file.write("\t".join(sline) + "\n")

            matched_length = str(len(matched_ids))
            ensembl_length = str(len(pathway_genes))

            print matched_length + " of " + ensembl_length + " IDs were found. "

            percent_matched = (float(len(matched_ids))/(float(len(pathway_genes))))*100
            print "percent matched: " + str(percent_matched)


        if float(matched_length) > (int(min_genes)-1 if min_genes else 4):
            used_pathways.append(pathway)

        else:
            print "Not enough matched genes...skipping eVIP for pathway"


    return pway_dict, used_pathways

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def run_eVIP_multiple_testing_pt1(infile=None, zscore_gct = None, out_directory=None, sig_info =None, c=None, r=None, num_reps=None,
         ie_filter=None,ie_col=None, i=None, allele_col=None, conn_null=None, conn_thresh=None,
         mut_wt_rep_rank_diff=None, use_c_pval=None, cell_id=None, plate_id=None, ref_allele_mode=None,
         x_thresh=None, y_thresh=None, annotate=None, by_gene_color=None, pdf=None, xmin=None,
         xmax=None, ymin=None, ymax=None, viz_ymin=None, viz_ymax=None, corr_val=None,sparkler_off=None,viz_off=None):

    #runs eVIP_corr and compare and returns the p values from compare
    #different sig_gctx for exp an z inputs used in viz

    sig_gctx_val = out_directory+ "/z_scores.gct"


    #run eVIP_corr.py
    # print('calculating correlations...')
    run_corr = eVIP_corr.run_main(input=infile,zscore_gct=zscore_gct, out_dir= out_directory)

    #run eVIP_compare.py
    # print('comparing...')
    run_compare = eVIP_compare.run_main(sig_info, out_directory+"/spearman_rank_matrix.gct", allele_col , out_directory+"/compare", r ,c , i , conn_null , ie_col, ie_filter, num_reps, cell_id , plate_id)

    #getting the p values from the pathway compare file
    mut_wt_rep_pvals_from_compare = []
    mut_wt_conn_null_pvals_from_compare = []
    wt_mut_rep_vs_wt_mut_conn_pvals_from_compare = []

    num_test_alleles = None

    with open(out_directory+"/compare.txt", "r") as compare_file:
        file_reader = csv.DictReader(compare_file, delimiter="\t")

        for row in file_reader:
            mut_wt_rep_pvals_from_compare.append(row['mut_wt_rep_pval'])
            mut_wt_conn_null_pvals_from_compare.append(row['mut_wt_conn_null_pval'])
            wt_mut_rep_vs_wt_mut_conn_pvals_from_compare.append(row['wt_mut_rep_vs_wt_mut_conn_pval'])


    #counting the number of mutations being tested
    with open(out_directory + "/compare.txt", "r") as compare_file:
        file_reader = csv.DictReader(compare_file, delimiter="\t")
        each_line = list(file_reader)
        num_test_alleles=len(each_line)

    return mut_wt_rep_pvals_from_compare, mut_wt_conn_null_pvals_from_compare, wt_mut_rep_vs_wt_mut_conn_pvals_from_compare,num_test_alleles

def run_eVIP_multiple_testing_pt2(infile=None, zscore_gct = None, out_directory=None, sig_info =None, c=None, r=None, num_reps=None,
         ie_filter=None,ie_col=None, i=None, allele_col=None, conn_null=None, conn_thresh=None,
         mut_wt_rep_rank_diff=None, use_c_pval=None, cell_id=None, plate_id=None, ref_allele_mode=None,
         x_thresh=None, y_thresh=None, annotate=None, by_gene_color=None, pdf=None, xmin=None,
         xmax=None, ymin=None, ymax=None, viz_ymin=None, viz_ymax=None, corr_val=None,mut_wt_rep_thresh=None,disting_thresh=None,sparkler_off=None,viz_off=None,cond_median_max_diff_thresh=None):

    sig_gctx_val = out_directory+ "/z_scores.gct"

    # print('predicting...')
    run_predict = eVIP_predict.run_main(out_directory+"/adj_compare.txt", out_directory+"/predict", conn_thresh,
                mut_wt_rep_thresh, mut_wt_rep_rank_diff,
                disting_thresh, use_c_pval,cond_median_max_diff_thresh)


    if not sparkler_off:
        # print "making sparkler plots..."
        run_sparkler = eVIP_sparkler.eVIP_run_main(out_directory+"/predict.txt", ref_allele_mode,
                y_thresh , x_thresh,
                use_c_pval,annotate, by_gene_color, pdf,
                xmin, xmax, ymin, ymax, out_directory+"/sparkler_plots")

    if not viz_off:
        if conn_null:
            null_conn = conn_null
        else:
            null_conn = out_directory + "/compare_conn_null.txt"

        """
        print "making visualizations..."
        run_viz = eVIP_viz.eVIP_run_main(pred_file= out_directory+"/predict.txt", sig_info = args.sig_info, gctx=out_directory+"/spearman_rank_matrix.gct",
                sig_gctx = sig_gctx_val, ref_allele_mode = args.ref_allele_mode, null_conn = null_conn,
                out_dir = out_directory+"/viz",ymin = args.viz_ymin, ymax= args.viz_ymax, allele_col = args.allele_col, use_c_pval = args.use_c_pval,
                 pdf = args.pdf, cell_id = args.cell_id, plate_id = args.plate_id, corr_val_str= args.corr_val)
        """
    # print "Pathway is DONE"

def run_eVIP(infile=None, zscore_gct = None, out_directory=None, sig_info =None, c=None, r=None, num_reps=None,
         ie_filter=None,ie_col=None, i=None, allele_col=None, conn_null=None, conn_thresh=None,
         mut_wt_rep_rank_diff=None, use_c_pval=None, cell_id=None, plate_id=None, ref_allele_mode=None,
         x_thresh=None, y_thresh=None, annotate=None, by_gene_color=None, pdf=None, xmin=None,
         xmax=None, ymin=None, ymax=None, viz_ymin=None, viz_ymax=None, corr_val=None,mut_wt_rep_thresh=None,disting_thresh=None,sparkler_off=None,viz_off=None,cond_median_max_diff_thresh=None):

    #different sig_gctx for exp an z inputs used in viz

    sig_gctx_val = out_directory+ "/z_scores.gct"
    # if args.zscore_gct :
    #     sig_gctx_val = args.zscore_gct


    # run eVIP_corr.py
    # print('calculating correlations...')
    run_corr = eVIP_corr.run_main(input=infile,zscore_gct=zscore_gct, out_dir= out_directory)

    # print('comparing...')
    run_compare = eVIP_compare.run_main(sig_info=sig_info, gctx = out_directory+"/spearman_rank_matrix.gct",
                allele_col = allele_col, o= out_directory+"/compare", r = r,
             c = c, i = i, conn_null = conn_null, ie_col = ie_col,
             ie_filter = ie_filter, num_reps = num_reps, cell_id = cell_id, plate_id = plate_id)

    # print('predicting...')
    run_predict = eVIP_predict.run_main(out_directory+"/compare.txt", out_directory+"/predict", conn_thresh,
                mut_wt_rep_thresh, mut_wt_rep_rank_diff,
                disting_thresh, use_c_pval,cond_median_max_diff_thresh)


    if not sparkler_off:
        # print "making sparkler plots..."
        run_sparkler = eVIP_sparkler.eVIP_run_main(pred_file = out_directory+"/predict.txt", ref_allele_mode=ref_allele_mode,
                y_thresh = y_thresh , x_thresh = x_thresh,
                use_c_pval= use_c_pval,annotate=annotate, by_gene_color= by_gene_color, pdf= pdf,
                xmin= xmin, xmax = xmax, ymin = ymin, ymax = ymax, out_dir = out_directory+"/sparkler_plots")

    if not viz_off:
        # print "making visualizations..."
        if conn_null:
            null_conn = conn_null
        else:
            null_conn = out_directory + "/compare_conn_null.txt"

        run_viz = eVIP_viz.eVIP_run_main(pred_file= out_directory+"/predict.txt", sig_info = sig_info, gctx=out_directory+"/spearman_rank_matrix.gct",
                sig_gctx = sig_gctx_val, ref_allele_mode = ref_allele_mode, null_conn = null_conn,
                out_dir = out_directory+"/viz",ymin = viz_ymin, ymax= viz_ymax, allele_col = allele_col, use_c_pval = use_c_pval,
                 pdf = pdf, cell_id = cell_id, plate_id = plate_id, corr_val_str= corr_val)

def padj(c_pval_type, num_test_alleles, pathway_list):
    # converting list of lists into one list
    c_pval_type = list(itertools.chain(*c_pval_type))
    #multiple testing corrections
    eVIPP_c_pval_type = robjects.r['p.adjust'](robjects.FloatVector(c_pval_type), "BH")

    new_cpval_list = str(eVIPP_c_pval_type).split()
    fixed_cpval_list = []

    for item in new_cpval_list:
        if item.startswith("["):
            pass
        else:
            fixed_cpval_list.append(item)

    #grouping pvalues from the same pathway together
    grouped_new_cpval_list = list(grouper(fixed_cpval_list, num_test_alleles, fillvalue=None))
    #zipping the pathway with the pvalues
    zipped = zip(grouped_new_cpval_list, pathway_list)

    return zipped

def grouper(iterable, n, fillvalue=None):
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)

def izip_longest(*args, **kwds):
    fillvalue = kwds.get('fillvalue')
    counter = [len(args) - 1]
    def sentinel():
        if not counter[0]:
            raise ZipExhausted
        counter[0] -= 1
        yield fillvalue
    fillers = repeat(fillvalue)
    iterators = [chain(it, sentinel(), fillers) for it in args]
    try:
        while iterators:
            yield tuple(map(next, iterators))
    except ZipExhausted:
        pass

def repeat(object, times=None):
    # repeat(10, 3) --> 10 10 10
    if times is None:
        while True:
            yield object
    else:
        for i in xrange(times):
            yield object

def chain(*iterables):
    # chain('ABC', 'DEF') --> A B C D E F
    for it in iterables:
        for element in it:
            yield element

class ZipExhausted(Exception):
    pass

def make_adj_compare_file(pathway_list, out_directory,eVIPP_adj_pways_mut_wt_rep_c_pvals_from_compare,eVIPP_adj_pways_mut_wt_conn_null_c_pvals_from_compare,eVIPP_adj_pways_wt_mut_rep_vs_wt_mut_conn_c_pvals_from_compare):

    column_headers = ["gene","mut","mut_rep", "wt_rep", "mut_wt_connectivity",
                      "wt", "cell_line", "mut_wt_rep_pval","mut_wt_conn_null_pval",
                      "wt_mut_rep_vs_wt_mut_conn_pval","kruskal_diff", "mut_wt_rep_c_pval",
                      "mut_wt_conn_null_c_pval","wt_mut_rep_vs_wt_mut_conn_c_pval"]

    for pathway in pathway_list:
        #new file
        adj_compare = open(out_directory + "/" + pathway + "_eVIPP_outputs/adj_compare.txt", "w")

        #writing header to new file
        file_writer = csv.DictWriter(adj_compare, delimiter="\t", fieldnames=column_headers)
        file_writer.writeheader()

        #opening the original compare file
        compare_file = open(out_directory + "/" +pathway+ "_eVIPP_outputs/compare.txt", "r")
        file_reader = csv.DictReader(compare_file, delimiter="\t")
        n = 0
        m = 0
        l = 0

        for line in file_reader:
            for item in eVIPP_adj_pways_mut_wt_rep_c_pvals_from_compare:
                if pathway in item:
                    line['mut_wt_rep_c_pval'] = item[0][n]
                    n+=1
            for item in eVIPP_adj_pways_mut_wt_conn_null_c_pvals_from_compare:
                if pathway in item:
                    line['mut_wt_conn_null_c_pval'] = item[0][m]
                    m += 1
            for item in eVIPP_adj_pways_wt_mut_rep_vs_wt_mut_conn_c_pvals_from_compare:
                if pathway in item:
                    line['wt_mut_rep_vs_wt_mut_conn_c_pval'] = item[0][l]
                    l += 1

            file_writer.writerow(line)
def summarize(pathway_list, output_file,out_dir):
    #creates a simple summary file of just variant predictions for each pathway

    mut_list = []
    pred_list = []
    #getting the list of tested mutations using the first pathway in the list
    with open(out_dir + "/" + pathway_list[0] + "_eVIPP_outputs/predict.txt") as first_pathway:
        for line in first_pathway:
            if line.startswith("gene"):
                continue
            else:
                splitLine = line.split()
                mut_list.append(str(splitLine[1]))

    header = "pathway" + "\t" + ("\t").join(mut_list) + "\n"
    output_file.write(header)

    for pathway in pathway_list:
        with open(out_dir + "/" + pathway + "_eVIPP_outputs/predict.txt") as this_predict:
            for line in this_predict:
                if line.startswith("gene"):
                    continue
                else:
                    splitLine = line.split()
                    pred_list.append(splitLine[14])

        output_file.write(pathway + "\t" + ("\t").join(pred_list)+"\n")
        #clearing predict list
        pred_list = []

def summarize_predict_files(pathway_list, output_file,out_dir):
    #combines predict files from each pathway

    # getting header using the first pathway in the list
    with open(out_dir + "/" + pathway_list[0] + "_eVIPP_outputs/predict.txt") as first_pathway:
        for line in first_pathway:
            if line.startswith("gene"):
                output_file.write("Pathway" +"\t"+ line)


    for pathway in pathway_list:
        with open(out_dir + "/" + pathway + "_eVIPP_outputs/predict.txt") as this_predict:
            for line in this_predict:
                if line.startswith("gene"):
                    continue
                else:
                    output_file.write(pathway +"\t"+ line)

    output_file.close()
