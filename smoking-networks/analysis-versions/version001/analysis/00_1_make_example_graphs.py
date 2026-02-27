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
out_dir = Path('results/00_1/')
out_dir.mkdir(exist_ok=True)


# \\\\
# \\\\
# Make Graphs
# \\\\
# \\\\


x = ["A", "B", "C", "D", "E", "F", "G","H", "I", "J", "K"]

edge_inds_A = [(0, 1), (0, 3), (1, 2), (3, 2),(3,4), (4, 6), (6, 5), (6, 4), (5, 4), (6,7), (7,8), (9,10)]
edges_A = [(x[a], x[b]) for a,b in edge_inds_A]

edge_inds_B = [(0, 1),(0,3), (1, 2), (3, 2), (4, 6), (6, 5), (6, 4), (5, 4), (6,7), (7,8), (9,10), (7,10)]
edges_B = [(x[a], x[b]) for a,b in edge_inds_B]

G_A = nx.Graph()
G_A.add_edges_from(edges_A)

G_B = nx.Graph()
G_B.add_edges_from(edges_B)


def plotG(G, G_other, title, ax):
    unique_edges_G = [(a,b) 
                    for a,b in G.edges()
                    if (b,a) not in G_other.edges() 
                        and (a,b) not in G_other.edges()
                    ]

    shared_edges_G = [(a,b) 
                    for a,b in G.edges()
                    if (b,a) in G_other.edges() 
                        or (a,b) in G_other.edges()
                    ]
    ax.axis("off")
    pos = nx.spring_layout(G, seed = 2026)
    nx.draw_networkx_nodes(G, ax = ax, pos = pos,node_size=1000, node_color = "#9c5ac2ff")
    nx.draw_networkx_labels(G, ax = ax, pos = pos, font_size=25,font_color = "white")
    nx.draw_networkx_edges(G, edgelist=shared_edges_G, edge_color="black", pos = pos, width = 3, ax = ax)
    nx.draw_networkx_edges(G, edgelist=unique_edges_G, edge_color= "orange", pos = pos, width = 6, ax = ax)
    ax.set_title(title)

fig, axes = plt.subplots(1,2,figsize = (14,7) )

plotG(G_A, G_B, title = None, ax = axes[0])
plotG(G_B, G_A, title = None, ax = axes[1])

# Add vertical line between subplots
fig.add_artist(plt.Line2D([0.5, 0.5], [0, 1], 
                          transform=fig.transFigure, 
                          color="grey", 
                          linestyle="--",
                          linewidth=2))
# plt.Line2D.set_linestyle(":")

# # Add manual legend
# from matplotlib.lines import Line2D
# legend_elements = [
#     Line2D([0], [0], color='orange', linewidth=4, label='Unique'),
#     Line2D([0], [0], color='black', linewidth=2, label='Shared')
# ]
# fig.legend(handles=legend_elements, loc='upper right', fontsize=12, frameon=True)
plt.tight_layout()
plt.savefig("results/00_1/shared-nodes-unq-edges.pdf")



