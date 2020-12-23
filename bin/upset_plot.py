########################################################################
# author: Alexis M. Thornton
#
# create a upset plot to see how much gene overlap there is between the
# different pathways from the pathway predictions
#
########################################################################
import pandas as pd
import json
import upsetplot
from upsetplot import from_contents
# pip install UpSetPlot==0.4.0
from matplotlib import pyplot as plt


def run(json_file,gmt,kallisto_spec,eVIPP_predict,output_name):

    spec_genes = pd.read_csv(kallisto_spec , sep="\t", index_col = "#gene_id").index.tolist()
    pathways = pd.read_csv(eVIPP_predict,sep="\t",index_col = "Pathway").index.tolist()

    if json_file:
        with open(json_file, "rb") as JSON:
            old_dict = json.load(JSON)
            gene_set_dict = {k: v for k, v in old_dict.iteritems() if v is not None}


    if gmt:
        with open(gmt) as gmt_:
            content = gmt_.read().splitlines()
            lines = [i.split("\t") for i in content]
            gene_set_dict = {item[0]: item[2:] for item in lines }
    #subset
    gene_set_dict = {k:v for (k,v) in gene_set_dict.items() if k in pathways}

    if len(gene_set_dict)>1:
        spec_dict = {}
        for k,v in gene_set_dict.items():
            spec_dict[k] = [i for i in v if i in spec_genes]

        e = from_contents(spec_dict)

        upsetplot.plot(e ,sort_by='cardinality', sort_categories_by='cardinality',show_counts=True )
        plt.savefig(output_name)
        plt.clf()


if __name__ == "__main__":
    main()
