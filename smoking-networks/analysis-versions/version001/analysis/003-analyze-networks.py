
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import networkx as nx
import pickle as pkl
import platform
import seaborn as sns

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

#\\
# look at edge weight distributions
#\\

# make pseudo-partial-correlation matrix
A = rcfgl[0][:,:,0]
B = rcfgl[0][:,:,1]

partial_A = abs(CovtoCor(A))
partial_B = abs(CovtoCor(B))

a, b = [np.log(x[x>0]) for x in [partial_A.ravel(), partial_B.ravel()]]

fig, ax = plt.subplots(figsize=(8, 4))
sns.kdeplot(a, ax=ax, color="steelblue", label="former",alpha = .5)
sns.kdeplot(b, ax=ax, color="orange", label="current", alpha = .5)
ax.legend()
plt.title("Distributions of log(partial-correlation) values where greater than zero")
plt.savefig('results/003/partial-corr-distr.pdf')

sum(partial_A.ravel()>np.exp(-6))
sum(partial_B.ravel()>np.exp(-6))

# make adjacency
truncation_threshold = .035 # .065 for 300
adj_all = MakeAdjMatrix_all(rcfgl,truncation_value= truncation_threshold, top_N = "all", names = "default")

# \\\\
# \\\\
# Look at node overlaps
# \\\\
# \\\\

def jaccard(s1, s2):
    """Intersection/union"""
    return  len(np.intersect1d(s1,s2))/len(np.union1d(s1,s2))

# metabolites kept after truncation/filtering
inds_kept = [x.columns.values for x in adj_all]

inds_intersect = np.intersect1d(*inds_kept)
# all 
inds_union = np.union1d(*inds_kept)
former  = np.isin(inds_union, inds_kept[0])
current = np.isin(inds_union, inds_kept[1])
    
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
upset_data = from_indicators(["former", "current"],pd.DataFrame({"former": former,"current": current,},index=inds_union))
up = UpSet(upset_data.copy(), subset_size="count", element_size=60)
fig = plt.figure(figsize = (20,20))
up.plot(fig = fig)
plt.title(f"Jaccard= {jaccard(*inds_kept).__round__(2)}\n for tau = {truncation_threshold}")
plt.savefig('results/003/former-current-metab-upset.pdf')
plt.close(fig)

# \\\\
# \\\\
# Plot Networks: 2x unique, and both
# \\\\
# \\\\
adj_former = adj_all[0].astype(int)
adj_current = adj_all[1].astype(int)
A_f = adj_former.loc[inds_intersect, inds_intersect].astype(bool)
A_c = adj_current.loc[inds_intersect, inds_intersect].astype(bool)

adj_both = A_c & A_f
adj_only_former =  ~A_c & A_f
adj_only_current =  A_c & ~A_f

def plot_graph(adj, rel_out_path, title, lab):
    # remove all unconnected nodes
    unconn = np.equal(adj.sum(1).values, 0) # rowsums
    adj = adj.iloc[~unconn, ~unconn]
    G = nx.from_pandas_adjacency(df=adj, create_using=nx.MultiGraph)
    nms = [rowDat[lab][n] for n in G.nodes()]
    labels = dict(zip(G.nodes(), nms))
    fig, ax = plt.subplots(figsize=(10, 10))
    pos = nx.spring_layout(G, seed=120349)
    nx.draw(G,pos,labels=labels,font_size=8,ax=ax)
    ax.set_title(title, pad=3)
    plt.savefig(rel_out_path)
    plt.close()

plot_graph(adj_both, 'results/003/both-graph.pdf', 'Both', lab = "sub_pathway")
plot_graph(adj_only_former, 'results/003/only-former-graph.pdf', 'Only in former smokers', lab = "sub_pathway")
plot_graph(adj_only_current, 'results/003/only-current-graph.pdf', 'Only in current smokers', lab = "sub_pathway")
plot_graph(adj_current, 'results/003/current-graph.pdf', 'Current smokers', lab = "sub_pathway")
plot_graph(adj_former, 'results/003/former-graph.pdf', 'Former smokers', lab = "sub_pathway")

plot_graph(adj_both, 'results/003/both-chem-name-graph.pdf', 'Both', lab = "chemical_name")
plot_graph(adj_only_former, 'results/003/only-former-chem-name-graph.pdf', 'Only in former smokers', lab = "chemical_name")
plot_graph(adj_only_current, 'results/003/only-current-chem-name-graph.pdf', 'Only in current smokers', lab = "chemical_name")
plot_graph(adj_current, 'results/003/current-chem-name-graph.pdf', 'Current smokers', lab = "chemical_name")
plot_graph(adj_former, 'results/003/former-chem-name-graph.pdf', 'Former smokers', lab = "chemical_name")


# \\\\
# \\\\
# Characterize main hubs
# \\\\
# \\\\

#\\\\
# current and former at argmax(degree(current))
#\\\\
adj = adj_current
hub = adj.sum(1).idxmax()
conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
inds = np.r_[conn_to_hub, hub]
hub_sub = adj.loc[inds, inds]

lab = "sub_pathway"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/current-hi-degree-sub-graph.pdf', f'Current smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)

lab = "chemical_name"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/current-hi-degree-sub-chem-name-graph.pdf', f'Current smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)


adj = adj_former
conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
inds = np.r_[conn_to_hub, hub]
hub_sub = adj.loc[inds, inds]

lab = "sub_pathway"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/former-at-current-hi-degree-sub-graph.pdf', f'Former smokers at current smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)

lab = "chemical_name"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/former-at-current-hi-degree-sub-chem-name-graph.pdf', f'Former smokers at current smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)



#\\\\
# current and former at argmax(degree(former))
#\\\\
adj = adj_former
hub = adj.sum(1).idxmax()
conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
inds = np.r_[conn_to_hub, hub]
hub_sub = adj.loc[inds, inds]

lab = "sub_pathway"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/former-hi-degree-sub-graph.pdf', f'Former smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)

lab = "chemical_name"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/former-hi-degree-sub-chem-name-graph.pdf', f'Former smokers at largest hub: {hub_metab}\nmaxdegree = {round(adj.loc[hub,:].sum())}', lab = lab)


adj = adj_current
conn_to_hub = adj.loc[hub, adj.loc[hub,:] >0 ].index
inds = np.r_[conn_to_hub, hub]
hub_sub = adj.loc[inds, inds]

lab = "sub_pathway"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/current-at-former-hi-degree-sub-graph.pdf', f'Current smokers at former smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)

lab = "chemical_name"
hub_metab = rowDat[lab].loc[hub] 
plot_graph(hub_sub, 'results/003/current-at-former-hi-degree-sub-chem-name-graph.pdf', f'Current smokers at former smokers largest hub: {hub_metab}\nmaxdegree={round(adj.loc[hub,:].sum())}', lab = lab)
