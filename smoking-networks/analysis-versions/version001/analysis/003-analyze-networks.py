
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

root= Path('/projects/canderson2@xsede.org/kechris-lab/smoking-networks')

try:
    os.chdir(root / 'analysis-versions/version001')
except OSError as e:
    print(f"Failed to set working directory: {e}", file=sys.stderr)
    sys.exit(1)

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

# functions for network analysis
from analysis.utils.myDstream_functions import *

import sklearn
iris = sklearn.datasets.load_iris(as_frame = True)['data']

cov = pd.DataFrame(np.cov(iris.T), index = iris.columns, columns = iris.columns)
prec = pd.DataFrame(np.linalg.inv(cov), index = iris.columns, columns = iris.columns)

partial_cor = np.abs(CovtoCor(prec))
adj, nms = MakeAdjMatrix(prec, .8, "all", 'default')
adj; nms

G = nx.from_numpy_array(adj, create_using=nx.MultiGraph)

fig = plt.figure(figsize=(10,10))
nx.draw_networkx(G)
plt.savefig('results/003/example-net.png')

# \\\\
# \\\\
# Visualize Network
# \\\\
# \\\\

G_former  = nx.from_pandas_adjacency(prec_former_f, create_using=nx.DiGraph)
G_current = nx.from_pandas_adjacency(prec_current_f, create_using=nx.DiGraph)

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
