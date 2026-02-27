### COMPONENT PAIRS ###
### COMPONENT PAIRS ###
### COMPONENT PAIRS ###
### COMPONENT PAIRS ###
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
out_dir = Path('results/005_2/')
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


# \\\\\
# \\\\\
# Find Component Pairs
# \\\\\
# \\\\\

# load components
from analysis.utils.myDstream_functions import MetabLoader

# load components
comp_curr = MetabLoader(path="results/005/curr-metabs/comp-metabs.txt").load(at=2, id=True, rowData=rowData)
comp_form = MetabLoader(path="results/005/form-metabs/comp-metabs.txt").load(at=2, id=True, rowData=rowData)

def jaccard(x,y):
    return len(np.intersect1d(x,y)) / len(np.union1d(x,y))    

res  = np.zeros((len(comp_curr), len(comp_form)))
for i in range(len(comp_curr)):
    for j in range(len(comp_form)):
        res[i,j] = jaccard(comp_curr[i], comp_form[j]) 

#\\\
# Create bipartite graph 
#\\\

G = nx.Graph()

# Add nodes for both sets
curr_nodes = [f"curr_{i+1}" for i in range(len(comp_curr))]
form_nodes = [f"form_{j+1}" for j in range(len(comp_form))]

G.add_nodes_from(curr_nodes, bipartite=0)
G.add_nodes_from(form_nodes, bipartite=1)

# Add edges with weights from the Jaccard matrix
for i in range(len(comp_curr)):
    for j in range(len(comp_form)):
        G.add_edge(curr_nodes[i], form_nodes[j], weight=res[i, j])

# Scale weights for better visualization
edges = G.edges()
weights = [G[u][v]['weight'] for u, v in edges]
edge_widths = [w*10 for w in weights]  

fig, ax = plt.subplots(figsize=(8, 10))
pos = nx.bipartite_layout(G, nodes=curr_nodes)
nx.draw_networkx_edges(G, pos, width=edge_widths, alpha=1, ax=ax)
nx.draw_networkx_nodes(G, pos, nodelist=curr_nodes, node_color='lightblue', node_size=800, label='Current', ax=ax)
nx.draw_networkx_nodes(G, pos, nodelist=form_nodes, node_color='lightcoral', node_size=800, label='Former', ax=ax)
nx.draw_networkx_labels(G, pos, font_size=10, ax=ax)
plt.legend()
plt.title("Current vs Former Components\nEdges weighted by Jaccard similarity")
plt.axis('off')
plt.tight_layout()
plt.savefig("results/005_2/component-pairs-graph.pdf")

