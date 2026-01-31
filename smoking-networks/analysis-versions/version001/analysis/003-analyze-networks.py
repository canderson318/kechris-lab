
import os
import re
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import networkx as nx
import pickle as pkl


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

# functions for network analysis
from analysis.utils.myDstream_functions import *

# set results directory
Path("results/003").mkdir(exist_ok=True)

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
# Process Network
# \\\\
# \\\\

truncation_threshold = .05

# make adjacency
adj_all = MakeAdjMatrix_all(rcfgl,truncation_value= truncation_threshold, top_N = "all", names = "default")

len(adj_all)
[x.shape for x in adj_all]

# \\\\
# \\\\
# Analyze Network
# \\\\
# \\\\

adj_former = adj_all[0].astype(int)
adj_current = adj_all[1].astype(int)

# metabolites kept after truncation/filtering
inds_kept = [x.columns.values for x in adj_all]

def jaccard(s1, s2):
    """Intersection/union"""
    return  len(np.intersect1d(s1,s2))/len(np.union1d(s1,s2))

jaccard(*inds_kept)

# metabolites in both
ind_intersect = np.intersect1d(*inds_kept)


both = adj_current.loc[ind_intersect, ind_intersect].astype(bool) & adj_former.loc[ind_intersect, ind_intersect].astype(bool)
only_former =  ~adj_current.loc[ind_intersect, ind_intersect].astype(bool) & adj_former.loc[ind_intersect, ind_intersect].astype(bool)
only_current =  adj_current.loc[ind_intersect, ind_intersect].astype(bool) & ~adj_former.loc[ind_intersect, ind_intersect].astype(bool)

only_former.sum(0).min()

def plot_graph(adj, rel_out_path, title):
    # remove all unconnected nodes
    unconn = np.equal(adj.sum(1).values, 0) # rowsums
    adj = adj.iloc[~unconn, ~unconn]
    G = nx.from_pandas_adjacency(df=adj, create_using=nx.MultiGraph)
    metabs = [rowDat.chemical_name[n] for n in G.nodes()]
    labels = dict(zip(G.nodes(), metabs))
    fig, ax = plt.subplots(figsize=(20, 20))
    pos = nx.spring_layout(G, seed=120349)
    nx.draw(G,pos,labels=labels,font_size=8,ax=ax)
    ax.set_title(title, pad=10)
    plt.savefig(rel_out_path)
    plt.close()

plot_graph(both, 'results/003/both-graph.pdf', 'Both')
plot_graph(only_former, 'results/003/only-former-graph.pdf', 'Only in former')
plot_graph(only_current, 'results/003/only-current-graph.pdf', 'Only in current')



print("\n•••••Former Smokers•••••")
NetworkPlotter(Adjacency_all, which = 1) 

print("\n•••••Current Smokers•••••")
NetworkPlotter(Adjacency_all, which = 2) 

PairNetworkPlotter(Adjacency_all, pair = [1, 2])

# \\\\
# \\\\
# xxx``
# \\\\
# \\\\
