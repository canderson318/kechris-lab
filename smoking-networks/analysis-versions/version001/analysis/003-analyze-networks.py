import os
import re
import sys
import numpy as np
import pandas as pd
import seaborn as sns
from pathlib import Path
import matplotlib.pyplot as plt
import networkx as nx

root= Path('/projects/canderson2@xsede.org/kechris-lab/smoking-networks')
RCFGL_path = root / 'RCFGL'
os.chdir(RCFGL_path)
try: 
    sys.path.remove("Python_functions")
except:
    print('Python_functions not in $PATH')
sys.path.insert(0, 'Python_functions')

from RCFGL import RFGL, RCFGL
from Dstream_functions import*

os.chdir(root / 'analysis-versions/version001')

# set results directory
Path("results/003").mkdir(exist_ok=True)

# \\\\
# \\\\
# Load Data 
# \\\\
# \\\\

rowDat = pd.read_csv('processed-data/001/rowData.csv', index_col = 0)

adj_dir = Path('results/002/RCFGL-output/adjacency/')
adj_files = os.listdir(adj_dir)
adj_names = [str(Path(x).stem) for x in adj_files]
adj = dict()
for nm, file in zip(adj_names, adj_files):
    adj[nm] = np.loadtxt(adj_dir/file, delimiter = ",")

prec_dir = Path('results/002/RCFGL-output/precision/')
prec_files = os.listdir(prec_dir)
prec_names = [str(Path(x).stem) for x in prec_files]
prec = dict()
for nm, file in zip(prec_names, prec_files):
    prec[nm] = np.loadtxt(prec_dir/file, delimiter = ",")

[x for x in adj.values()][0].shape

# \\\\
# \\\\
# Visualize Network
# \\\\
# \\\\

# G = nx.DiGraph()
# G.add_nodes_from(rowDat.metab_id)
prec_current, prec_former = [x for x in prec.values()]

def process(mat, thresh):
    """Format data for graphing with thresholded edges"""
    mat_f = mat.copy()
    mat_f[mat_f < thresh] = 0
    for i,j in (mat.shape[0], mat.shape[1]):
        mat[i,j] = 0
    df = pd.DataFrame(
        mat_f,
        index=rowDat.metab_id,
        columns=rowDat.metab_id,
    )
    # remove zero-weight edges
    df = df.loc[:, (df != 0).any(axis=0)]
    df = df.loc[(df != 0).any(axis=1), :]
    return df

thr = .65
G_former = nx.from_pandas_adjacency(process(prec_former, thresh  = thr), create_using=nx.DiGraph, )
G_current = nx.from_pandas_adjacency(process(prec_current, thresh  = thr), create_using=nx.DiGraph)

fig = plt.figure(figsize=(30, 30))
pos = nx.spring_layout(G_former, seed=120349)
nx.draw(G_former, pos, with_labels=True, font_size=8)
plt.title("Former Smokers")
plt.savefig("results/003/former-prec-graph.pdf")
plt.close()


print("\n•••••Never Smokers•••••")
NetworkPlotter(Adjacency_all, which = 1) 

print("\n•••••Former Smokers•••••")
NetworkPlotter(Adjacency_all, which = 2) 

print("\n•••••Current Smokers•••••")
NetworkPlotter(Adjacency_all, which = 3) 

PairNetworkPlotter(Adjacency_all, pair = [1, 2])
# \\\\
# \\\\
# xxx``
# \\\\
# \\\\
