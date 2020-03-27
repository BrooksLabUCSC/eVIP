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


def run(json_file,kallisto_spec,eVIPP_predict,output_name):

    spec_genes = pd.read_csv(kallisto_spec , sep="\t", index_col = "#gene_id").index.tolist()
    pathways = pd.read_csv(eVIPP_predict,sep="\t",index_col = "Pathway").index.tolist()

    with open(json_file) as f:
        gene_set_dict = json.load(f)

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
