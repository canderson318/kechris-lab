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
import networkx as nx
import pickle as pkl
import platform
import seaborn as sns
import warnings

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
Path("results/005").mkdir(exist_ok=True)

# \\\\
# \\\\
# Load Data 
# \\\\
# \\\\

rowDat = pd.read_csv('processed-data/001/rowData.csv', index_col = 0)

rcfgl = pkl.load(open('results/003/RCFGL-output/RCFGL.pkl','rb'))

len(rcfgl)
prec_array, _, _ = rcfgl
prec_array.shape

# make precisions dataframes
nms = ['curr', "form"]
prec_dict = {nms[i]: pd.DataFrame(prec_array[:,:,i], index = rowDat.metab_id, columns = rowDat.metab_id, dtype = float) for i in range(prec_array.shape[2])}

# \\\\\
# \\\\\
#  Make graphs and adjacency
# \\\\\
# \\\\\

from analysis.utils.myDstream_functions import Adjacency, MakeAdjMatrix, CovtoCor

def adjGraphs(prec_list, truncation_value, top_N):
    adjacencies = [
        Adjacency(theta=prec_list[i], truncation_value=truncation_value, top_N=top_N)
        for i in range(len(prec_list))
    ]
    graphs = [
        nx.from_pandas_adjacency(df=adjacencies[i], create_using=nx.MultiGraph) 
        for i in range(len(prec_list))
    ]
    if not np.array_equal(np.array(graphs[0].nodes()), adjacencies[0].index.values):
        raise ValueError("Graph nodes do not equal adjacency index/columns")
    return graphs[0], adjacencies[0],   graphs[1], adjacencies[1]

# make graphs for all edges > 0
tau = truncation_value = 0.04
top_N  = 'all'
theta = list(prec_dict.values())[0]

G_curr, adj_curr,  G_form, adj_form = adjGraphs(list(prec_dict.values()), tau, top_N)

adj_curr.to_csv('results/005/curr-adj.csv')
adj_form.to_csv('results/005/form-adj.csv')

# \\\\
# \\\\
# List hubs
# \\\\
# \\\\
def save_top_hubs(G, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    conn_comp = list(nx.connected_components(G))
    lngs = [len(x) for x in conn_comp]
    top_hub_ids = [conn_comp[i] for i in np.argsort(lngs)[-5:]]
    top_hub_chems = [list(rowDat.chemical_name[rowDat.metab_id.isin(ids)].values) for ids in top_hub_ids]
    top_hub_super = [list(rowDat.super_pathway[rowDat.metab_id.isin(ids)].values) for ids in top_hub_ids]
    top_hub_sub = [list(rowDat.sub_pathway[rowDat.metab_id.isin(ids)].values) for ids in top_hub_ids]
    sep = '•'
    with open(out_dir/'metabs.txt', 'w', encoding="UTF8") as f:
        for super, sub, chem in zip(top_hub_super,top_hub_sub, top_hub_chems):
            pathways = [f"{sp}{sep}{sb}{sep}{ch}" for sp, sb, ch in zip(super, sub, chem)]
            print("||".join(pathways),'\n', file=f)



save_top_hubs(G = G_form, out_dir = 'results/005/form-hub-metabs/')
save_top_hubs(G = G_curr, out_dir = 'results/005/curr-hub-metabs/')


# \\\\
# \\\\
# Look at node overlaps
# \\\\
# \\\\


def jaccard(s1, s2):
    """Intersection/union"""
    return  len(np.intersect1d(s1,s2))/len(np.union1d(s1,s2))

# metabolites kept after truncation/filtering
inds_kept = [x.columns.values for x in (adj_form, adj_curr)]

inds_intersect = np.intersect1d(*inds_kept)

# all 
inds_union = np.union1d(*inds_kept)
form  = np.isin(inds_union, inds_kept[0])
curr = np.isin(inds_union, inds_kept[1])
     
from upsetplot import UpSet, from_indicators
upset_data = from_indicators(["form", "curr"],pd.DataFrame({"form": form,"curr": curr,},index=inds_union))
up = UpSet(upset_data.copy(), subset_size="count", element_size=60)
fig = plt.figure(figsize = (20,20))
up.plot(fig = fig)
plt.title(f"Jaccard= {jaccard(*inds_kept).__round__(2)}\n for tau = {tau}\nand/or top N = {top_N}")
plt.savefig('results/005/form-curr-metab-upset.pdf')
plt.close(fig)




# # \\\\
# # \\\\
# # Plot Networks: 2x unique, and both
# # \\\\
# # \\\\

# inds_kept = [x.columns.values for x in (adj_form, adj_curr)]

# inds_intersect = np.intersect1d(*inds_kept)

# # all 
# inds_union = np.union1d(*inds_kept)

# A_f = adj_form.loc[inds_intersect, inds_intersect].astype(bool)
# A_c = adj_curr.loc[inds_intersect, inds_intersect].astype(bool)

# adj_both = A_c & A_f
# adj_only_form =  ~A_c & A_f
# adj_only_curr =  A_c & ~A_f

# def plot_graph(adj, rel_out_path, title, lab, size = (10,10)):
#     # remove all unconnected nodes
#     unconn = np.equal(adj.sum(1).values, 0) # rowsums
#     adj = adj.iloc[~unconn, ~unconn]
#     G = nx.from_pandas_adjacency(df=adj, create_using=nx.MultiGraph)
#     nms = rowDat[rowDat.metab_id.isin(G.nodes())][lab] 
#     labels = dict(zip(G.nodes(), nms))
#     fig, ax = plt.subplots(figsize=size)
#     # pos = nx.spring_layout(G, seed=120349)
#     pos = nx.circular_layout(G)
#     nx.draw(G,pos,labels=labels,font_size=8,ax=ax)
#     ax.set_title(title, pad=3)
#     plt.savefig(rel_out_path)
#     plt.close()


# plot_graph(adj_both, 'results/005/both-graph.pdf', 'Both', lab = "sub_pathway", size = (20,20))
# plot_graph(adj_only_form, 'results/005/only-form-graph.pdf', 'Only in form smokers', lab = "sub_pathway")
# plot_graph(adj_only_curr, 'results/005/only-curr-graph.pdf', 'Only in curr smokers', lab = "sub_pathway")
# plot_graph(adj_curr, 'results/005/curr-graph.pdf', 'curr smokers', lab = "sub_pathway")
# plot_graph(adj_form, 'results/005/form-graph.pdf', 'form smokers', lab = "sub_pathway")

# plot_graph(adj_both, 'results/005/both-chem-name-graph.pdf', 'Both', lab = "chemical_name", size = (20,20))
# plot_graph(adj_only_form, 'results/005/only-form-chem-name-graph.pdf', 'Only in form smokers', lab = "chemical_name")
# plot_graph(adj_only_curr, 'results/005/only-curr-chem-name-graph.pdf', 'Only in curr smokers', lab = "chemical_name")
# plot_graph(adj_curr, 'results/005/curr-chem-name-graph.pdf', 'curr smokers', lab = "chemical_name")
# plot_graph(adj_form, 'results/005/form-chem-name-graph.pdf', 'form smokers', lab = "chemical_name")


# # \\\\
# # \\\\
# # Characterize main hubs
# # \\\\
# # \\\\

# #\\\\
# # curr and form at argmax(degree(curr))
# #\\\\
# adj = adj_curr
# hub = adj.sum(1).idxmax()
# conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index.values
# inds = np.concatenate([conn_to_hub, [hub]])
# hub_sub = adj.loc[inds, inds]

# lab = "sub_pathway"
# hub_metab = rowDat[rowDat.metab_id==hub][lab].item()
# plot_graph(hub_sub, 'results/005/curr-hi-degree-sub-graph.pdf', f'curr smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)

# lab = "chemical_name"
# plot_graph(hub_sub, 'results/005/curr-hi-degree-sub-chem-name-graph.pdf', f'curr smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)


# adj = adj_form
# conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
# inds = np.concatenate([conn_to_hub, [hub]])
# hub_sub = adj.loc[inds, inds]

# lab = "sub_pathway"
# plot_graph(hub_sub, 'results/005/form-at-curr-hi-degree-sub-graph.pdf', f'form smokers at curr smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)

# lab = "chemical_name"
# plot_graph(hub_sub, 'results/005/form-at-curr-hi-degree-sub-chem-name-graph.pdf', f'form smokers at curr smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)



# #\\\\
# # curr and form at argmax(degree(form))
# #\\\\
# adj = adj_form
# hub = adj.sum(1).idxmax()
# conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
# inds = np.concatenate([conn_to_hub, [hub]])
# hub_sub = adj.loc[inds, inds]

# lab = "sub_pathway"
# hub_metab = rowDat[rowDat.metab_id==hub][lab].item()
# plot_graph(hub_sub, 'results/005/form-hi-degree-sub-graph.pdf', f'form smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)

# lab = "chemical_name"
# plot_graph(hub_sub, 'results/005/form-hi-degree-sub-chem-name-graph.pdf', f'form smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)


# adj = adj_curr
# conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
# inds = np.concatenate([conn_to_hub, [hub]])
# hub_sub = adj.loc[inds, inds]

# lab = "sub_pathway"
# plot_graph(hub_sub, 'results/005/curr-at-form-hi-degree-sub-graph.pdf', f'curr smokers at form smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)

# lab = "chemical_name"
# plot_graph(hub_sub, 'results/005/curr-at-form-hi-degree-sub-chem-name-graph.pdf', f'curr smokers at form smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)
