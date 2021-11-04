import dash
import dash_core_components as dcc
import dash_html_components as html
import dash_table

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import plotly.express as px

import pandas as pd
import numpy as np

import os

import sys

inp_dir = None

#os.chdir('../')

# Give the user options to use the raw p-value (by default), or to use the corrected p-value
use_c_pval = False
try:
    if len(sys.argv) >= 2 and sys.argv[1] == '--use_c_pval':
        use_c_pval = True
        if len(sys.argv) >= 3 and sys.argv[2] == '-inp-dir':
            inp_dir = sys.argv[3]
    elif len(sys.argv) >= 2 and sys.argv[1] == '-inp-dir':
        inp_dir = sys.argv[2]
    else:
        inp_dir = os.path.join(os.getcwd(), r'tutorial_files/eVIP2_out') # Default
except:
    pass

# Some basic styling features
external_stylesheets = ['https://codepen.io/chriddyp/pen/bWLwgP.css']

# This creates the web app which will be rendered
app = dash.Dash(__name__, external_stylesheets=external_stylesheets, suppress_callback_exceptions = True)

# The navbar will contain buttons to navigate between each page
navbar_items = ["Overall Impact", "Pathway-specific Impact", "Help"]

# This loads in the gene data used in the overall impact section into a pandas dataframe 
predict_path = None
if inp_dir == None:
    predict_path = os.path.join(os.getcwd(), r'tutorial_files/eVIP2_out/eVIP_out/predict.txt')
else:
    predict_path = os.path.join(inp_dir, r'eVIP_out/predict.txt')
    
df = pd.read_csv(predict_path, delimiter='\t')

if not use_c_pval: 
    df["mut_wt_rep_c_pval"] = df["mut_wt_rep_pval"]

for col in df.columns:
    if col not in ['gene','mut', 'mut_rep', 'wt_rep', 'wt_mut_rep_vs_wt_mut_conn_c_pval', "mut_wt_rep_c_pval",  'prediction']:
        df = df.drop(col, axis=1)

# Extract some features which will be used when plotting graphs in the overall impact section
df["impact_pval"] = np.log10(df["wt_mut_rep_vs_wt_mut_conn_c_pval"]) * -1
yvals = []
for dat in zip(df["mut_rep"].values, df["wt_rep"].values, df["mut_wt_rep_c_pval"].values):
    mut_rep = dat[0]
    wt_rep = dat[1]
    mut_wt_rep_c_pval = dat[2]
    if mut_rep > wt_rep:
        yvals.append(-np.log(mut_wt_rep_c_pval))
    else:
        yvals.append(np.log(mut_wt_rep_c_pval))
df["impact_direction"] = yvals

# This will read in the data for each gene variant
pathway_dict = dict()
mut_path_dict = dict()

for gene_var in set(df["mut"]):
    #path = os.path.join(os.getcwd(), 'testing-data')
    path = os.getcwd()
    if inp_dir != None:
        path = inp_dir
    else: 
        path = os.path.join(path, r'tutorial_files/eVIP2_out')
    
    path = os.path.join(path, 'eVIPP_out') 
    path = os.path.join(path, gene_var)
    wt_path = os.path.join(path, 'wt_specific')
    mut_path = os.path.join(path, 'mutation_specific')

    # Get the wild type and mutant data path
    wt_path = os.path.join(wt_path, 'eVIPP_combined_predict_files.txt')
    mut_path = os.path.join(mut_path, 'eVIPP_combined_predict_files.txt')

    # Try reading the mutation data - if it doesn't exist, then it will be empty
    try: 
        path_df = pd.read_csv(wt_path, delimiter="\t")

        if not use_c_pval: 
            path_df["mut_wt_rep_c_pval"] = path_df["mut_wt_rep_pval"]
            path_df["wt_mut_rep_vs_wt_mut_conn_c_pval"] = path_df["wt_mut_rep_vs_wt_mut_conn_pval"]
        
        for col in path_df.columns:
            if col not in ["Pathway", "wt_mut_rep_vs_wt_mut_conn_c_pval", "prediction", "mut_wt_rep_c_pval", "mut_rep", "wt_rep"]:
                path_df.drop(col, inplace=True, axis=1)
        pathway_dict[gene_var] = path_df
        pathway_dict[gene_var]["impact_pval"] = np.log10(pathway_dict[gene_var]["wt_mut_rep_vs_wt_mut_conn_c_pval"]) * -1
        
        yvals = []
        for dat in zip(path_df["mut_rep"].values, path_df["wt_rep"].values, path_df["mut_wt_rep_c_pval"].values):
            mut_rep = dat[0]
            wt_rep = dat[1]
            mut_wt_rep_c_pval = dat[2]
            if mut_rep > wt_rep:
                yvals.append(-np.log(mut_wt_rep_c_pval)) # todo: add option to use raw pval
            else:
                yvals.append(np.log(mut_wt_rep_c_pval))
        pathway_dict[gene_var]["impact_direction"] = yvals

    except Exception as e:
        pathway_dict[gene_var] = None
    
    
    try: 
        mut_path_df = pd.read_csv(mut_path, delimiter="\t")

        if not use_c_pval: 
            mut_path_df["mut_wt_rep_c_pval"] = mut_path_df["mut_wt_rep_pval"]
            mut_path_df["wt_mut_rep_vs_wt_mut_conn_c_pval"] = mut_path_df["wt_mut_rep_vs_wt_mut_conn_pval"]

        for col in mut_path_df.columns:
            if col not in ["Pathway", "wt_mut_rep_vs_wt_mut_conn_c_pval", "prediction", "mut_wt_rep_c_pval", "mut_rep", "wt_rep"]:
                mut_path_df.drop(col, inplace=True, axis=1)
        mut_path_dict[gene_var] = mut_path_df
        mut_path_dict[gene_var]["impact_pval"] = np.log10(mut_path_dict[gene_var]["wt_mut_rep_vs_wt_mut_conn_c_pval"]) * -1
        
        yvals = []
        for dat in zip(mut_path_df["mut_rep"].values, mut_path_df["wt_rep"].values, mut_path_df["mut_wt_rep_c_pval"].values):
            mut_rep = dat[0]
            wt_rep = dat[1]
            mut_wt_rep_c_pval = dat[2]
            if mut_rep > wt_rep:
                yvals.append(-np.log(mut_wt_rep_c_pval))
            else:
                yvals.append(np.log(mut_wt_rep_c_pval))
        mut_path_dict[gene_var]["impact_direction"] = yvals

    except Exception as e:
        mut_path_dict[gene_var] = None


for gene in set(df["gene"]):
    specific_gene_df = df[df["gene"] == gene]
    x_data = specific_gene_df["impact_pval"].values
    y_data = specific_gene_df["impact_direction"].values
    z_data = specific_gene_df["prediction"].values # This will contain the hover text
    zz_data = specific_gene_df['mut'].values 

    fig = go.Figure()
    #fig.add_scatter(x=x_data, y=y_data, mode='markers', hoverinfo='text', hovertext=df['mut'].values, showlegend=False)

    for tup in zip(x_data, y_data, z_data, zz_data):
        # Solution to 0, 0: do the line thing
        z = tup[2]
        zz = tup[3]
        hover_data = ["", zz]
        if z == "GOF":
            fig.add_trace(go.Scatter(x=[0, tup[0]], y=[0, tup[1]], hovertext=hover_data, hoverinfo='text', marker=dict(size=17, opacity=[0, 1]), marker_color='rgba(202, 0, 32, 0.6)',
                        line_color='rgba(202, 0, 32, 0.6)', hoverlabel=dict(bgcolor='black'), showlegend=False))
        elif z == "LOF":
            fig.add_trace(go.Scatter(x=[0, tup[0]], y=[0, tup[1]], hovertext=hover_data, hoverinfo='text', marker=dict(size=17, opacity=[0, 1]), marker_color='rgba(5, 113, 176, 0.6)', 
                        line_color='rgba(5, 113, 176, 0.6)', hoverlabel=dict(bgcolor='black'), showlegend=False))
        elif z == "COF":
            fig.add_trace(go.Scatter(x=[0, tup[0]], y=[0, tup[1]], hovertext=hover_data, hoverinfo='text', marker=dict(size=17, opacity=[0, 1]), marker_color='rgba(106, 61, 154, 0.6)',
                        line_color='rgba(106, 61, 154, 0.6)', hoverlabel=dict(bgcolor='black'), showlegend=False))
        else:
            fig.add_trace(go.Scatter(x=[0, tup[0]], y=[0, tup[1]], marker=dict(size=17, opacity=[0, 1]), hovertext=hover_data, hoverinfo='text', marker_color='rgba(0, 0, 0, 0.6)', line_color='rgba(0, 0, 0, 0.6)', hoverlabel=dict(bgcolor='black'), showlegend=False))

    # increase circle size, change alpha value
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                marker=dict(size=10, color='#CA0020'),
                legendgroup='GOF', showlegend=True, name='GOF'))

    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
            marker=dict(size=10, color='#0571B0'),
            legendgroup='LOF', showlegend=True, name='LOF'))

    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
            marker=dict(size=10, color='#6A3D9A'),
            legendgroup='COF', showlegend=True, name='COF'))
    
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
            marker=dict(size=10, color='#000000'),
            legendgroup='Neutral', showlegend=True, name='Neutral'))

    fig.update_layout(height=600, width=800, xaxis_title="-log10 (corrected p-val)", yaxis_title="Impact Direction Score", showlegend=True, yaxis_range=[-4, 4], xaxis_range=[0, 4])

    if not use_c_pval:
        fig.update_layout(xaxis_title="-log10 (pval)")
    
    fig.write_html(os.path.join(inp_dir, str(gene) + '.html'))
    
    
    def convert_overall_column(col_name):
        if col_name == 'mut':
            return 'Mutant'
        elif col_name == 'wt_mut_rep_vs_wt_mut_conn_c_pval':
            return 'Impact p-val'
        return 'Prediction'
    
    rendered_df = pd.concat([specific_gene_df['mut'], specific_gene_df['wt_mut_rep_vs_wt_mut_conn_c_pval'], specific_gene_df['prediction']], axis=1)
    rendered_df = rendered_df.sort_values(by=['wt_mut_rep_vs_wt_mut_conn_c_pval'])
    rendered_table = dash_table.DataTable(id='table', columns=[{"name" : convert_overall_column(i), "id": i} for i in rendered_df.columns], 
                                data=rendered_df.to_dict('records'))



for gene_var in set(df["mut"]):
    pathway_table = None # wild type specific table 
    pathway_graph = None
    mut_table = None # mutation specific table
    mut_graph = None

    def convert_overall_column(col_name):
        if col_name == 'Pathway':
            return 'Pathway'
        elif col_name == 'wt_mut_rep_vs_wt_mut_conn_c_pval':
            return 'Impact p-val'
        return 'Prediction'

    try: 
        rendered_df = pathway_dict[gene_var]

        x_data = rendered_df["impact_pval"].values
        y_data = rendered_df["impact_direction"].values
        z_data = rendered_df["prediction"].values
        zz_data = rendered_df['Pathway'].values 

        fig = go.Figure()
        #fig.add_scatter(x=x_data, y=y_data, mode='markers', hoverinfo='text', hovertext=rendered_df['Pathway'].values, showlegend=False)

        for tup in zip(x_data, y_data, z_data, zz_data):
            z = tup[2]
            zz = tup[3]

            hover_data = ["", zz]

            if z == "GOF":
                fig.add_trace(go.Scatter(x=[0, tup[0]], y=[0, tup[1]], marker=dict(size=17, opacity=[0, 1]), marker_color='rgba(202, 0, 32, 0.6)', hovertext=hover_data, hoverinfo='text', hoverlabel=dict(bgcolor='black'), line_color='rgba(202, 0, 32, 0.6)', showlegend=False))
            elif z == "LOF":
                fig.add_trace(go.Scatter(x=[0, tup[0]], y=[0, tup[1]], marker=dict(size=17, opacity=[0, 1]), marker_color='rgba(5, 113, 176, 0.6)', hovertext=hover_data, hoverinfo='text', line_color='rgba(5, 113, 176, 0.6)', hoverlabel=dict(bgcolor='black'), showlegend=False))
            elif z == "COF":
                fig.add_trace(go.Scatter(x=[0, tup[0]], y=[0, tup[1]], marker=dict(size=17, opacity=[0, 1]), marker_color='rgba(106, 61, 154, 0.6)', hovertext=hover_data, hoverinfo='text', line_color='rgba(106, 61, 154, 0.6)', hoverlabel=dict(bgcolor='black'), showlegend=False))
            else:
                fig.add_trace(go.Scatter(x=[0, tup[0]], y=[0, tup[1]], marker=dict(size=17, opacity=[0, 1]), marker_color='rgba(0, 0, 0, 0.6)', hovertext=hover_data, hoverinfo='text', line_color='rgba(0, 0, 0, 0.6)', hoverlabel=dict(bgcolor='black'), showlegend=False))

        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                marker=dict(size=10, color='#CA0020'),
                legendgroup='GOF', showlegend=True, name='GOF'))

        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                marker=dict(size=10, color='#0571B0'),
                legendgroup='LOF', showlegend=True, name='LOF'))

        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                marker=dict(size=10, color='#6A3D9A'),
                legendgroup='COF', showlegend=True, name='COF'))
        
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                marker=dict(size=10, color='#000000'),
                legendgroup='Neutral', showlegend=True, name='Neutral'))

        fig.update_layout(height=600, width=800, xaxis_title="-log10 (corrected p-val)", yaxis_title="Impact Direction Score", showlegend=True, 
                            yaxis_range=[-4, 4], xaxis_range=[0, 4])

        if not use_c_pval:
            fig.update_layout(xaxis_title="-log10 (pval)")

        pathway_graph = html.Div(dcc.Graph(id='pathway_graph', figure=fig))

        rendered_df = rendered_df.sort_values(by=['wt_mut_rep_vs_wt_mut_conn_c_pval'])
        rendered_df = pd.concat([rendered_df['Pathway'], rendered_df['wt_mut_rep_vs_wt_mut_conn_c_pval'], rendered_df['prediction']], axis=1)
        pathway_table = dash_table.DataTable(id='pathway-table', columns=[{"name" : convert_overall_column(i), "id": i} for i in rendered_df.columns], 
                        data=rendered_df.to_dict('records'))
    except Exception as e:
        pathway_table = None
        pathway_graph = None
        pass

    if pathway_graph != None:
        fig.write_html(os.path.join(inp_dir, str(gene_var) + '_wt_specific.html'))

    try:
        rendered_df = mut_path_dict[gene_var] # fix mut_path_dict


        x_data = rendered_df["impact_pval"].values
        y_data = rendered_df["impact_direction"].values
        z_data = rendered_df["prediction"].values
        zz_data = rendered_df['Pathway'].values 
        
        fig = go.Figure()
        # fig.add_scatter(x=x_data, y=y_data, mode='markers', hoverinfo='text', hovertext=rendered_df['Pathway'].values, showlegend=False)

        for tup in zip(x_data, y_data, z_data, zz_data):
            zz = tup[3]

            hover_data = ["", zz]

            z = tup[2]
            if z == "GOF":
                fig.add_trace(go.Scatter(x=[0, tup[0]], y=[0, tup[1]], marker=dict(size=17, opacity=[0, 1]), marker_color='rgba(202, 0, 32, 0.6)', hovertext=hover_data, hoverinfo='text', line_color='rgba(202, 0, 32, 0.6)', hoverlabel=dict(bgcolor='black'), showlegend=False) )
            elif z == "LOF":
                fig.add_trace(go.Scatter(x=[0, tup[0]], y=[0, tup[1]], marker=dict(size=17, opacity=[0, 1]), marker_color='rgba(5, 113, 176, 0.6)', hovertext=hover_data, hoverinfo='text', line_color='rgba(5, 113, 176, 0.6)', hoverlabel=dict(bgcolor='black'), showlegend=False))
            elif z == "COF":
                fig.add_trace(go.Scatter(x=[0, tup[0]], y=[0, tup[1]], marker=dict(size=17, opacity=[0, 1]), marker_color='rgba(106, 61, 154, 0.6)', hovertext=hover_data, hoverinfo='text', line_color='rgba(106, 61, 154, 0.6)', hoverlabel=dict(bgcolor='black'), showlegend=False))
            else:
                fig.add_trace(go.Scatter(x=[0, tup[0]], y=[0, tup[1]], marker=dict(size=17, opacity=[0, 1]), marker_color='rgba(0, 0, 0, 0.6)', hovertext=hover_data, hoverinfo='text', line_color='rgba(0, 0, 0, 0.6)', hoverlabel=dict(bgcolor='black'), showlegend=False))


        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                marker=dict(size=10, color='#CA0020'),
                legendgroup='GOF', showlegend=True, name='GOF'))

        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                marker=dict(size=10, color='#0571B0'),
                legendgroup='LOF', showlegend=True, name='LOF'))

        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                marker=dict(size=10, color='#6A3D9A'),
                legendgroup='COF', showlegend=True, name='COF'))
        
        fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers',
                marker=dict(size=10, color='#000000'),
                legendgroup='Neutral', showlegend=True, name='Neutral'))

        fig.update_layout(height=600, width=800, xaxis_title="-log10 (corrected p-val)", yaxis_title="Impact Direction Score", showlegend=True, 
                        yaxis_range=[-4, 4], xaxis_range=[0, 4])

        if not use_c_pval:
            fig.update_layout(xaxis_title="-log10 (pval)")


        mut_graph = html.Div(dcc.Graph(id='mut_graph', figure=fig))

        rendered_df = rendered_df.sort_values(by=['wt_mut_rep_vs_wt_mut_conn_c_pval'])
        rendered_df = pd.concat([rendered_df['Pathway'], rendered_df['wt_mut_rep_vs_wt_mut_conn_c_pval'], rendered_df['prediction']], axis=1)
        mut_table = dash_table.DataTable(id='mut_table', columns=[{"name": convert_overall_column(i), "id": i} for i in rendered_df.columns], 
                data=rendered_df.to_dict('records'))
    except Exception as e:
        mut_table = None
        mut_graph = None
        pass

    if mut_graph != None:
        fig.write_html(os.path.join(inp_dir, str(gene_var) + '_mutation_specific.html'))
        '''
        wt_h2 = html.H6("WT-specific Sparkler Plot")
        mut_h2 = html.H6("Mutation-specific Sparkler Plot")

        wt_h2_tab = html.H6("WT-specific Table")
        mut_h2_tab = html.H6("Mutation-specific Table")

        gene = df[df["mut"] == gene_var]["gene"]
        gene_h2 = html.H4(gene)

        # FIND NUMBER OF DIFERENTIALLY EXPRESSED GENES
        path = os.path.join(os.getcwd(), 'testing-data')
        path = os.path.join(path, 'kallisto_files') 
        mut_path = os.path.join(path, 'combined_kallisto_abundance_genes_filtered_transformed_' + gene_var + '_mutspec.tsv')
        wt_path = os.path.join(path, 'combined_kallisto_abundance_genes_filtered_transformed_transformed_' + gene_var + '_wtspec.tsv')
        
        mut_diff = None
        num_diff_mut = None
        wt_diff = None
        num_wt_mut = None

        try: 
            mut_diff = pd.read_csv(mut_path, delimiter='\t')
            num_diff_mut = len(mut_diff.index)
        except Exception as e:
            num_diff_mut = None

        try: 
            wt_diff = pd.read_csv(wt_path, delimiter='\t')
            num_wt_mut = len(wt_diff.index)
        except Exception as e:
            num_diff_mut = None

        # EVEN IF THE # OF DIFFERENTIALLY EXPRESSED GENES EXIST, THEY WILL ONLY BE DISPLAYED IF THE 
        # SUMMARY IS IN THE WT_SPECIFIC OR MUTATION SPECIFIC FOLDERS

        num_diff_wt_render = html.H6(str(num_wt_mut) + " Wildtype-specific differentially expressed genes")
        num_diff_mut_render = html.H6(str(num_diff_mut) + " Mutation-specific differentially expressed genes")

        
        prediction = df[df['mut'] == gene_var]['prediction']
        prediction = str(prediction.values[0])

        gene_var_header = html.H6(gene_var + ": " + prediction)

        
        if mut_table == None:
            return html.Div(
                children=[
                    gene_h2, gene_var_header, html.Hr(), num_diff_wt_render, wt_h2, pathway_graph, wt_h2_tab, pathway_table, num_diff_mut_render
                ]
            )
        elif pathway_table == None:
            return html.Div(
                children=[
                    gene_h2, gene_var_header, html.Hr(), num_diff_wt_render, num_diff_mut_render, mut_h2, mut_graph, mut_h2_tab, mut_table
                ]
            )
        

        return html.Div(
            children=[
                gene_h2, gene_var_header, html.Hr(), num_diff_wt_render, wt_h2, pathway_graph, wt_h2_tab, pathway_table, num_diff_mut_render, mut_h2, mut_graph, mut_h2_tab, mut_table
            ]
        )'''
        