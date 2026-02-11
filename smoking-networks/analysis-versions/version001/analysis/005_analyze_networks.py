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

rcfgl = pkl.load(open('results/002/RCFGL-output/RCFGL.pkl','rb'))

len(rcfgl)
prec_array, _, _ = rcfgl
prec_array.shape



# \\\\
# \\\\
# Plot Networks: 2x unique, and both
# \\\\
# \\\\
A_f = adj_form.loc[inds_intersect, inds_intersect].astype(bool)
A_c = adj_curr.loc[inds_intersect, inds_intersect].astype(bool)

adj_both = A_c & A_f
adj_only_form =  ~A_c & A_f
adj_only_curr =  A_c & ~A_f

def plot_graph(adj, rel_out_path, title, lab, size = (10,10)):
    # remove all unconnected nodes
    unconn = np.equal(adj.sum(1).values, 0) # rowsums
    adj = adj.iloc[~unconn, ~unconn]
    G = nx.from_pandas_adjacency(df=adj, create_using=nx.MultiGraph)
    nms = [rowDat[lab][n] for n in G.nodes()]
    labels = dict(zip(G.nodes(), nms))
    fig, ax = plt.subplots(figsize=size)
    # pos = nx.spring_layout(G, seed=120349)
    pos = nx.circular_layout(G)
    nx.draw(G,pos,labels=labels,font_size=8,ax=ax)
    ax.set_title(title, pad=3)
    plt.savefig(rel_out_path)
    plt.close()

plot_graph(adj_both, 'results/005/both-graph.pdf', 'Both', lab = "sub_pathway", size = (20,20))
plot_graph(adj_only_form, 'results/005/only-form-graph.pdf', 'Only in form smokers', lab = "sub_pathway")
plot_graph(adj_only_curr, 'results/005/only-curr-graph.pdf', 'Only in curr smokers', lab = "sub_pathway")
plot_graph(adj_curr, 'results/005/curr-graph.pdf', 'curr smokers', lab = "sub_pathway")
plot_graph(adj_form, 'results/005/form-graph.pdf', 'form smokers', lab = "sub_pathway")

plot_graph(adj_both, 'results/005/both-chem-name-graph.pdf', 'Both', lab = "chemical_name", size = (20,20))
plot_graph(adj_only_form, 'results/005/only-form-chem-name-graph.pdf', 'Only in form smokers', lab = "chemical_name")
plot_graph(adj_only_curr, 'results/005/only-curr-chem-name-graph.pdf', 'Only in curr smokers', lab = "chemical_name")
plot_graph(adj_curr, 'results/005/curr-chem-name-graph.pdf', 'curr smokers', lab = "chemical_name")
plot_graph(adj_form, 'results/005/form-chem-name-graph.pdf', 'form smokers', lab = "chemical_name")


# \\\\
# \\\\
# Characterize main hubs
# \\\\
# \\\\

#\\\\
# curr and form at argmax(degree(curr))
#\\\\
adj = adj_curr
hub = adj.sum(1).idxmax()
conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
inds = np.r_[conn_to_hub, hub]
hub_sub = adj.loc[inds, inds]

lab = "sub_pathway"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/005/curr-hi-degree-sub-graph.pdf', f'curr smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)

lab = "chemical_name"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/005/curr-hi-degree-sub-chem-name-graph.pdf', f'curr smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)


adj = adj_form
conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
inds = np.r_[conn_to_hub, hub]
hub_sub = adj.loc[inds, inds]

lab = "sub_pathway"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/005/form-at-curr-hi-degree-sub-graph.pdf', f'form smokers at curr smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)

lab = "chemical_name"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/005/form-at-curr-hi-degree-sub-chem-name-graph.pdf', f'form smokers at curr smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)



#\\\\
# curr and form at argmax(degree(form))
#\\\\
adj = adj_form
hub = adj.sum(1).idxmax()
conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
inds = np.r_[conn_to_hub, hub]
hub_sub = adj.loc[inds, inds]

lab = "sub_pathway"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/005/form-hi-degree-sub-graph.pdf', f'form smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)

lab = "chemical_name"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/005/form-hi-degree-sub-chem-name-graph.pdf', f'form smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)


adj = adj_curr
conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
inds = np.r_[conn_to_hub, hub]
hub_sub = adj.loc[inds, inds]

lab = "sub_pathway"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/005/curr-at-form-hi-degree-sub-graph.pdf', f'curr smokers at form smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)

lab = "chemical_name"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/005/curr-at-form-hi-degree-sub-chem-name-graph.pdf', f'curr smokers at form smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)
