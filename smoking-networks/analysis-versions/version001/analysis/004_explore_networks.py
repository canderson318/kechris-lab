### EXPLORE NETWORKS ###
### EXPLORE NETWORKS ###
### EXPLORE NETWORKS ###
### EXPLORE NETWORKS ###

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
Path("results/004").mkdir(exist_ok=True)

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

#\\
## make adjacency
#\\

from analysis.utils.myDstream_functions import Adjacency, pseudoPartialCorr, MakeAdjMatrix, CovtoCor

truncation_value = .035
top_N  = 'all'

adj_form = Adjacency(theta = prec_array[:,:,1], truncation_value=truncation_value, top_N = top_N)
adj_form=adj_form.astype(int)

# # check my function is correct
# adj_formx, _ = MakeAdjMatrix(theta = prec_array[:,:,0], truncation_value=truncation_value, top_N = top_N, names = 'none')
# np.array_equal( adj_form.to_numpy(), adj_formx)

adj_curr = Adjacency(prec_array[:,:,0],truncation_value=truncation_value, top_N = top_N)
adj_curr=adj_curr.astype(int)

G_form = nx.from_pandas_adjacency(df=adj_form, create_using=nx.MultiGraph) # Undirected Graphs
G_curr = nx.from_pandas_adjacency(df=adj_curr, create_using=nx.MultiGraph)

# \\\\
# \\\\
# Calculate network metrics
# \\\\
# \\\\

#\\
# look at edge weight distributions
#\\

partial_form = abs(pseudoPartialCorr(prec_array[:,:,1]))
partial_curr= abs(pseudoPartialCorr(prec_array[:,:,0]))

# transform to less skewed for plotting
form, curr = (partial_form.ravel(), partial_curr.ravel())

fig, ax = plt.subplots(figsize=(8, 4))
sns.kdeplot(form, ax=ax, color="steelblue", label="form",alpha = .5,clip = (None, 1) )
sns.kdeplot(curr, ax=ax, color="orange", label="curr", alpha = .5, clip = (None, 1) )
ax.legend()
plt.title("Distributions of partial-correlation")
plt.savefig('results/004/partial-corr-distr.pdf')

form, curr = [np.log(x[x>0]) for x in [partial_form.ravel(), partial_curr.ravel()]]
fig, ax = plt.subplots(figsize=(8, 4))
sns.kdeplot(form, ax=ax, color="steelblue", label="form",alpha = .9,clip = (None, 0) )
sns.kdeplot(curr, ax=ax, color="orange", label="curr", alpha = .9, clip = (None, 0) )
ax.legend()
plt.title("Distributions of log(partial-correlation>0)")
plt.savefig('results/004/log-partial-corr-distr.pdf')


#\\
# make table of 
#:: connectivity, edge density, sparsity, topoligy, count nodes, count edges,...
#\\

def network_metrix(G):
    adj = nx.adjacency_matrix(G).toarray().squeeze()
    n_node = adj.shape[0]
    n_edge = (adj>0).sum()
    n_sparse = (adj==0).sum()
    dens_edge = n_edge/np.multiply(*adj.shape)
    dens_sparse = n_sparse/np.multiply(*adj.shape)
    # connectivity
    conn_comp = list(nx.connected_components(G)) # islands
    largest_conn_comp = max([len(x) for x in conn_comp]) # largest island
    smallest_conn_comp = min([len(x) for x in conn_comp]) # smallest
    num_isol = len(list(nx.isolates(G))) # unconnected nodes
    degs = np.array([x for x in dict(G.degree()).values()])
    deg_mean = degs.mean()
    deg_max = degs.max()
    # num connected components
    # output
    return { "N Nodes": n_node,"N Edges": n_edge,"N Sparse": n_sparse, "Average Degree": deg_mean.round(3),"Max Degree": deg_max.round(3), "Edge Density": dens_edge.round(3),"Sparse Density":dens_sparse.round(3), "Largest Connected Component": largest_conn_comp,  "Smallest Connected Component": smallest_conn_comp,  "Number Connected Components": conn_comp.__len__(),  "Number Isolated Nodes": num_isol, }

metric_form, metric_curr = (network_metrix(G_form), network_metrix(G_curr))
df = pd.DataFrame({'': metric_form.keys(), 'Former': metric_form.values(),'Current': metric_curr.values()})

# save to pdf
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman', 'DejaVu Serif']

fig, ax = plt.subplots(figsize=(10, 5))  # Adjust size as needed
ax.axis('tight')
ax.axis('off')
table = ax.table(cellText=df.values, colLabels=df.columns, 
                 loc='center', cellLoc='right')
table.auto_set_font_size(False)
table.set_fontsize(12)
plt.savefig('results/004/network-metrics.pdf', bbox_inches='tight', pad_inches=0.1)

# L = nx.normalized_laplacian_matrix(G_form)
# e = np.linalg.eigvals(L.toarray()).ravel()
# sum(e<1e-5)
# sns.histplot(e, bins=30);
# plt.close()

#\\\
# degree densities
#\\\
deg_form = [x for x in dict(G_form.degree()).values()]
deg_curr = [x for x in dict(G_curr.degree()).values()]

df = pd.DataFrame({
    "degree": np.concatenate([deg_form, deg_curr]),
    "group": (["form"] * len(deg_form)) + (["curr"] * len(deg_curr)),
})

fig, ax = plt.subplots(figsize=(8, 5))
sns.histplot(data=df,x="degree",hue="group",fill=True,alpha=0.5,common_norm=False,multiple="dodge",  ax=ax,kde=True, )
ax.set_title("Graph Degrees")
ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
plt.savefig("results/004/degree-densities.pdf")

#\\\
# shortest path densities
#\\\
def short_paths(G):
    for s in G.nodes():
        lengths = nx.single_source_shortest_path_length(G, s)
        yield from lengths.values()

short_paths_form = [x for x in short_paths(G_form)]
short_paths_curr = [x for x in short_paths(G_curr)]

df = pd.DataFrame({
    "length": np.concatenate([short_paths_form, short_paths_curr]), 
    'group': (["form"] * len(short_paths_form)) + (["curr"] * len(short_paths_curr))
})

fig, ax = plt.subplots(figsize=(8, 5))
sns.histplot(data=df,x="length",hue="group",fill=True,alpha=0.5,common_norm=False,multiple="dodge",  ax=ax, kde=True,kde_kws={"bw_adjust":2})
ax.set_title("Shortest Path Lengths")
ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
plt.savefig("results/004/shortest-path-densities.pdf")


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
    top_hub_inds = [list(conn_comp[i]) for i in np.sort(lngs)[-5:]]
    
    top_hub_chems = [list(rowDat.chemical_name[inds].values) for inds in top_hub_inds]
    top_hub_super = [list(rowDat.super_pathway[inds].values) for inds in top_hub_inds]
    top_hub_sub = [list(rowDat.sub_pathway[inds].values) for inds in top_hub_inds]
    sep = '•'
    with open(out_dir/'metabs.txt', 'w', encoding="UTF8") as f:
        for super, sub, chem in zip(top_hub_super,top_hub_sub, top_hub_chems):
            pathways = [f"{sp}{sep}{sb}{sep}{ch}" for sp, sb, ch in zip(super, sub, chem)]
            print("||".join(pathways),'\n', file=f)

save_top_hubs(G_form, 'results/004/form-hub-metabs/')
save_top_hubs(G_curr, 'results/004/curr-hub-metabs/')


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
    
# fig, ax1 = plt.subplots()
# # Left y-axis: Jaccards (0–1)
# ax1.plot(x[:len(jaccards)], jaccards)
# ax1.set_ylim(0, 1)
# ax1.set_ylabel("Jaccard")
# # Right y-axis: Shapes (0–800)
# ax2 = ax1.twinx()
# ax2.plot(x[:len(shapes)], shapes, c = 'orange')
# ax2.set_ylim(0, 800)
# ax2.set_ylabel("Shapes")
# ax1.set_xlabel("x")
# plt.show()
# plt.close(fig)


from upsetplot import UpSet, from_indicators
upset_data = from_indicators(["form", "curr"],pd.DataFrame({"form": form,"curr": curr,},index=inds_union))
up = UpSet(upset_data.copy(), subset_size="count", element_size=60)
fig = plt.figure(figsize = (20,20))
up.plot(fig = fig)
plt.title(f"Jaccard= {jaccard(*inds_kept).__round__(2)}\n for tau = {truncation_value}\nand/or top N = {top_N}")
plt.savefig('results/004/form-curr-metab-upset.pdf')
plt.close(fig)




