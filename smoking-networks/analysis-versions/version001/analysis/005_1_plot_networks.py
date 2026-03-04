### ANALYZE NETWORKS ###
### ANALYZE NETWORKS ###
### ANALYZE NETWORKS ###
### ANALYZE NETWORKS ###
import os
import sys
import numpy as np
import subprocess as sp
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx
import pickle as pkl
import platform
import seaborn as sns
import warnings
from pprint import pprint

warnings.filterwarnings('ignore', message='.*chained assignment.*')

# assign highest level directory
op_sys = platform.system()

if op_sys == 'Linux':
    root=Path.home() / '/projects/canderson2@xsede.org/kechris-lab/smoking-networks/'
elif op_sys == "Darwin": 
    root=Path.home() / 'Documents/school/local-kechris-lab/kechris-lab/smoking-networks'
else:
    raise SystemError("System Not Identifiable")

try:
    os.chdir(root / 'analysis-versions/version001')
except OSError as e:
    print(f"Failed to set working directory: {e}", file=sys.stderr)
    sys.exit(1)

# set results directory
out_dir = Path('results/005_1/')
out_dir.mkdir(exist_ok=True)

# \\\\
# \\\\
# Load Data 
# \\\\
# \\\\

rowData = pd.read_csv('processed-data/001/rowData.csv', index_col = 0)

#  load graphs and adjacencies
def load_adjGraphs(paths: list):
    """Load adjacencies from last script and make graphs"""
    adjacencies = [pd.read_csv(p, index_col= 0) for p in paths]
    graphs = [
        nx.from_pandas_adjacency(df=adjacencies[i], create_using=nx.MultiGraph) 
        for i in range(len(adjacencies))
    ]
    if not np.array_equal(np.array(graphs[0].nodes()), adjacencies[0].index.values):
        raise ValueError("Graph nodes do not equal adjacency index/columns")
    return graphs[0], adjacencies[0],   graphs[1], adjacencies[1]

paths = ['results/005/curr-adj.csv','results/005/form-adj.csv']
G_curr, adj_curr,  G_form, adj_form = load_adjGraphs(paths)


# load components

# extract_metabs  <-  function(list){
#     lapply(list, function(S){
#         M = unlist(str_split(S, "\\|\\|"))
#         m = lapply(M, function(s){
#             str_split(s, '•') %>% 
#                 lapply(., `[`, 3)  %>% # grab third element, chemical name
#                 lapply(., function(x){
#                     str_trim(x)
#                 })
#         }) %>% 
#         unlist()
#     })
# }


def load_metabs(path, at, id = True):
    def _get_metab_id(comp_list):
        out = [rowData.metab_id.values[rowData.chemical_name.isin(component)] for component in comp_list]
        return out
    lines = open(path, 'r').readlines()
    metabs = [  [y.split("•")[at] for y in x.split("||")] for x in lines ]
    if id:
        return _get_metab_id(metabs)
    else: 
        return metabs

from analysis.utils.myDstream_functions import MetabLoader

# load components
comp_curr = MetabLoader(path="results/005/curr-metabs/comp-metabs.txt").load(at=2, id=True, rowData=rowData)
comp_form = MetabLoader(path="results/005/form-metabs/comp-metabs.txt").load(at=2, id=True, rowData=rowData)

# \\\\
# \\\\
# Plot Networks: 2x unique, and both
# \\\\
# \\\\


def plot_graph(adj, other_adj, out_dir, nm, title, size=(15,15), fntsz = 12, filter_to_common = True):
    # remove all unconnected nodes
    unconn = np.equal(adj.sum(1).values, 0)
    adj = adj.iloc[~unconn, ~unconn]
    
    # Get the intersection of nodes that exist in both matrices
    common_nodes = adj.index.intersection(other_adj.index)
    
    # alphabetize nodes
    common_nodes = pd.Index.sort_values(common_nodes)
    
    # save unique nodes
    unique_nodes = adj.index.difference(other_adj.index).to_numpy()

    if filter_to_common:
        # Filter both to common nodes
        adj = adj.loc[common_nodes, common_nodes]
        other_adj = other_adj.loc[common_nodes, common_nodes]

    # make graphs
    G = nx.from_pandas_adjacency(df=adj, create_using=nx.Graph)
    G_other = nx.from_pandas_adjacency(df=other_adj, create_using=nx.Graph)
    
    # unique to G edges
    unique_edges_list = [(u, v) for u, v in G.edges() if (u, v) not in G_other.edges() and (v, u) not in G_other.edges()]
    # shared edges
    shared_edges_list = [(u, v) for u, v in G.edges() if (u, v) in G_other.edges() or (v, u) in G_other.edges()]
    
    G_unique_nodes = [n for n in unique_nodes if n in G.nodes()]

    for lab in ["sub_pathway", "chemical_name"]:
        path = Path(out_dir)/f"{lab}-{nm}"
        
        unique_nms = rowData[rowData.metab_id.isin(G_unique_nodes)][lab]
        unique_labels = dict(zip(G_unique_nodes, unique_nms))
        
        all_nms = rowData[rowData.metab_id.isin(G.nodes())][lab]
        all_labels = dict(zip(G.nodes(), all_nms))

        
        fig, ax = plt.subplots(figsize=size)
        ax.axis("off")
        pos = nx.circular_layout(G)
        
        # draw G_other as base graph
        nx.draw_networkx_nodes(G, pos, ax=ax,  node_color="#b979deff", edgecolors = "white",node_size=500 )
        nx.draw_networkx_labels(G, pos, ax=ax, labels = all_labels, font_size=fntsz , font_color="black")
        nx.draw_networkx_labels(G, pos, ax=ax, labels=unique_labels, font_size=fntsz , font_color="orange")
        
        # overlay edges unique to G_other in orange
        nx.draw_networkx_edges(G, pos, edgelist=shared_edges_list, edge_color="black", width=2, ax=ax)
        nx.draw_networkx_edges(G, pos, edgelist=unique_edges_list, edge_color="orange", width=2, ax=ax)
        
        legend_elements = [
            Line2D([0], [0], color="orange", linewidth=2, label=f"unique"),
            Line2D([0], [0], color="black", linewidth=2, label=f"shared")
            ]
        ax.legend(handles=legend_elements, loc="best")
        ax.set_title(title, pad=3)
        plt.savefig(path)
        plt.close()

# metabolites kept after truncation/filtering
inds_curr, inds_form = [x.columns.values for x in (adj_curr, adj_form)]

# shared metabolites
inds_intersect = np.intersect1d(*(inds_curr, inds_form))

# adjacencies at shared nodes
adj_shared_curr = adj_curr.loc[inds_intersect, inds_intersect].astype(bool)
adj_shared_form = adj_form.loc[inds_intersect, inds_intersect].astype(bool)

#\\\ 
#\\\ Shared Nodes
#\\\ 
plot_graph(adj = adj_shared_curr,other_adj=adj_shared_form, out_dir = out_dir, nm =  'node-shared-curr-graph.pdf', title =  'Current Smokers At Shared Nodes',  size = (20,20))
plot_graph(adj = adj_shared_form, other_adj =  adj_shared_curr, out_dir =out_dir,nm = 'node-shared-form-graph.pdf', title = 'Former Smokers At Shared Nodes',  size = (20,20))

# (adj_shared_form & ~adj_shared_curr).sum().sum() / 2


#\\\ 
#\\\ Shared Edges
#\\\ 

# Edge same/distinctness for shared nodes
adj_edge_both = adj_shared_curr & adj_shared_form
# edges unique to current
adj_edge_distinct_curr =  adj_shared_curr & ~adj_shared_form
# edges unique to former
adj_edge_distinct_form =  ~adj_shared_curr & adj_shared_form

# plot edge level sameness/distinctness
plot_graph(adj_edge_both,adj_edge_both, out_dir, 'edge-shared-graph.pdf', 'Edges in Both',  size = (20,20))
plot_graph(adj_edge_distinct_curr,adj_edge_distinct_curr,  out_dir,'edge-only-curr-graph.pdf', 'Edges Only in Current Smokers')
plot_graph(adj_edge_distinct_form,adj_edge_distinct_form, out_dir,'edge-only-form-graph.pdf', 'Edges Only in Former Smokers')

#\\\ 
#\\\ Distinct Nodes
#\\\ 

# nodes unique to current
inds_distinct_curr = np.setdiff1d(inds_curr, inds_form)
# nodes unique to former
inds_distinct_form= np.setdiff1d(inds_form, inds_curr)

# adjacencies at distinct nodes
adj_node_distinct_curr = adj_curr.loc[inds_distinct_curr, inds_distinct_curr].astype(bool)
adj_node_distinct_form = adj_form.loc[inds_distinct_form, inds_distinct_form].astype(bool)

# plot edge level sameness/distinctness
# adj=adj_node_distinct_curr
# out_dir=out_dir
# nm='node-only-curr-graph.pdf'
# title='Nodes Only in Current Smokers'

plot_graph(adj=adj_node_distinct_curr,other_adj = adj_node_distinct_form, out_dir=out_dir,nm='node-only-curr-graph.pdf', title='Nodes Only in Current Smokers')
plot_graph(adj = adj_node_distinct_form,other_adj = adj_node_distinct_curr,out_dir =  out_dir,nm = 'node-only-form-graph.pdf',title =  'Nodes Only in Former Smokers')


#\\\
# plot components
#\\\
def jaccard(s1, s2):
    """Intersection/union"""
    return  len(np.intersect1d(s1,s2))/len(np.union1d(s1,s2))


# plot params
label_sizes = [13]*4 + [13]*2 + [12]*4
plot_sizes = [(10,8)]*4 + [(12,10)]*2 + [(17,15)]*4


out_dir = Path("results/005_1/component-graphs/current/")
out_dir.mkdir(exist_ok=True, parents=True)
for i  in range(len(comp_curr)):
    _adj_curr = adj_curr.loc[comp_curr[i], comp_curr[i]]
    if i in [2,5]: # use adj at components for these because they are so similar
    # if jaccard(comp_curr[i],comp_form[i]) > .8: # use adj at components for these because they are so similar
        _other_adj = adj_form.loc[comp_form[i],comp_form[i]]
        filter_to_common = False
    else:
        _other_adj = adj_form
        filter_to_common = True
    plot_graph(adj = _adj_curr, other_adj = _other_adj, out_dir = out_dir,nm =  f"{i}.pdf", 
               title = f"Component {i+1}", 
               fntsz = label_sizes[i], size = plot_sizes[i],
               filter_to_common = filter_to_common)

out_dir = Path("results/005_1/component-graphs/former/")
out_dir.mkdir(exist_ok=True)
for i in range(len(comp_form)):
    _adj_form = adj_form.loc[comp_form[i], comp_form[i]]
    if i in [2,5]: # use adj at components for these because they are so similar
    # if jaccard(comp_curr[i],comp_form[i]) > .9 : # use adj at components for these because they are so similar
        _other_adj = adj_curr.loc[comp_curr[i],comp_curr[i]]
        filter_to_common = False
    else:
        _other_adj = adj_curr
        filter_to_common = True
    plot_graph(adj = _adj_form,other_adj = _other_adj, out_dir= out_dir,nm =  f"{i}.pdf", 
               title = f"Component {i+1}", 
               fntsz = label_sizes[i],size=plot_sizes[i],
               filter_to_common = filter_to_common)

# i = 2
# lab = "chemical_name"
# size = (10,10)

# a_comp_curr = comp_curr[i]
# _adj_curr = adj_curr.loc[a_comp_curr, a_comp_curr]
# adj = _adj_curr
# other_adj = adj_form
# out_dir = out_dir
# nm =  f"{i}.pdf"
# title = f"Component {i+1}"

# a_comp_form = comp_form[i]
# _adj_form = adj_form.loc[a_comp_form, a_comp_form]
# adj = _adj_form
# other_adj = adj_curr
# out_dir= out_dir
# nm =  f"{i}.pdf"
# title = f"Component {i+1}"
