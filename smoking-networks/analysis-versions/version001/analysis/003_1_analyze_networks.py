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

#\\
## make adjacency
#\\

from analysis.utils.myDstream_functions import Adjacency, pseudoPartialCorr, MakeAdjMatrix, CovtoCor

truncation_value = .035
top_N  = 'all'
adj_form = Adjacency(theta = prec_array[:,:,0], truncation_value=truncation_value, top_N = top_N)
adj_form=adj_form.astype(int)

# # check my function is correct
# adj_formx, _ = MakeAdjMatrix(theta = prec_array[:,:,0], truncation_value=truncation_value, top_N = top_N, names = 'none')
# np.array_equal( adj_form.to_numpy(), adj_formx)

adj_curr = Adjacency(prec_array[:,:,1],truncation_value=truncation_value, top_N = top_N)
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

partial_form = abs(pseudoPartialCorr(prec_array[:,:,0]))
partial_curr= abs(pseudoPartialCorr(prec_array[:,:,1]))

# transform to less skewed for plotting
form, curr = (partial_form.ravel(), partial_curr.ravel())

fig, ax = plt.subplots(figsize=(8, 4))
sns.kdeplot(form, ax=ax, color="steelblue", label="form",alpha = .5)
sns.kdeplot(curr, ax=ax, color="orange", label="curr", alpha = .5)
ax.legend()
plt.title("Distributions of partial-correlation")
plt.savefig('results/003/partial-corr-distr.pdf')

form, curr = [np.log(x[x>0]) for x in [partial_form.ravel(), partial_curr.ravel()]]
fig, ax = plt.subplots(figsize=(8, 4))
sns.kdeplot(form, ax=ax, color="steelblue", label="form",alpha = .5)
sns.kdeplot(curr, ax=ax, color="orange", label="curr", alpha = .5)
ax.legend()
plt.title("Distributions of log(partial-correlation>0)")
plt.savefig('results/003/log-partial-corr-distr.pdf')


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
    # Average path length: average distance between any two nodes in a graph
    #+++
    out = {
        'N Nodes': n_node,'N Edges': n_edge,'N Sparse': n_sparse,
        'Average Degree': deg_mean.round(3),'Max Degree': deg_max.round(3),
        'Edge Density': dens_edge.round(3),'Sparse Density':dens_sparse.round(3),
        'Largest Connected Component': largest_conn_comp, 
        'Smallest Connected Component': smallest_conn_comp, 
        'Number Isolated Components': num_isol,
        }
    return out

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
plt.savefig('results/003/network-metrics.pdf', bbox_inches='tight', pad_inches=0.1)



# degree densities
deg_form = [x for x in dict(G_form.degree()).values()]
deg_curr = [x for x in dict(G_curr.degree()).values()]

df = pd.DataFrame({
    "degree": np.concatenate([deg_form, deg_curr]),
    "group": (["form"] * len(deg_form)) + (["curr"] * len(deg_curr)),
})

fig, ax = plt.subplots(figsize=(8, 5))
sns.histplot(data=df,x="degree",hue="group",fill=True,alpha=0.5,common_norm=False,multiple="dodge",  ax=ax,)
ax.set_title("Graph Degrees")
ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
plt.savefig("results/003/degree-densities.pdf")


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
plt.savefig('results/003/form-curr-metab-upset.pdf')
plt.close(fig)

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

plot_graph(adj_both, 'results/003/both-graph.pdf', 'Both', lab = "sub_pathway", size = (20,20))
plot_graph(adj_only_form, 'results/003/only-form-graph.pdf', 'Only in form smokers', lab = "sub_pathway")
plot_graph(adj_only_curr, 'results/003/only-curr-graph.pdf', 'Only in curr smokers', lab = "sub_pathway")
plot_graph(adj_curr, 'results/003/curr-graph.pdf', 'curr smokers', lab = "sub_pathway")
plot_graph(adj_form, 'results/003/form-graph.pdf', 'form smokers', lab = "sub_pathway")

plot_graph(adj_both, 'results/003/both-chem-name-graph.pdf', 'Both', lab = "chemical_name", size = (20,20))
plot_graph(adj_only_form, 'results/003/only-form-chem-name-graph.pdf', 'Only in form smokers', lab = "chemical_name")
plot_graph(adj_only_curr, 'results/003/only-curr-chem-name-graph.pdf', 'Only in curr smokers', lab = "chemical_name")
plot_graph(adj_curr, 'results/003/curr-chem-name-graph.pdf', 'curr smokers', lab = "chemical_name")
plot_graph(adj_form, 'results/003/form-chem-name-graph.pdf', 'form smokers', lab = "chemical_name")


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
plot_graph(hub_sub, 'results/003/curr-hi-degree-sub-graph.pdf', f'curr smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)

lab = "chemical_name"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/curr-hi-degree-sub-chem-name-graph.pdf', f'curr smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)


adj = adj_form
conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
inds = np.r_[conn_to_hub, hub]
hub_sub = adj.loc[inds, inds]

lab = "sub_pathway"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/form-at-curr-hi-degree-sub-graph.pdf', f'form smokers at curr smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)

lab = "chemical_name"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/form-at-curr-hi-degree-sub-chem-name-graph.pdf', f'form smokers at curr smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)



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
plot_graph(hub_sub, 'results/003/form-hi-degree-sub-graph.pdf', f'form smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)

lab = "chemical_name"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/form-hi-degree-sub-chem-name-graph.pdf', f'form smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)


adj = adj_curr
conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
inds = np.r_[conn_to_hub, hub]
hub_sub = adj.loc[inds, inds]

lab = "sub_pathway"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/curr-at-form-hi-degree-sub-graph.pdf', f'curr smokers at form smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)

lab = "chemical_name"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/curr-at-form-hi-degree-sub-chem-name-graph.pdf', f'curr smokers at form smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)
