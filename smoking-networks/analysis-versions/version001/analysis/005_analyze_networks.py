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

rowData = pd.read_csv('processed-data/001/rowData.csv', index_col = 0)

rcfgl = pkl.load(open('results/003/RCFGL-output/RCFGL.pkl','rb'))

len(rcfgl)
prec_array, _, _ = rcfgl
prec_array.shape

# make precisions dataframes
nms = ['curr', "form"]
prec_dict = {nms[i]: pd.DataFrame(prec_array[:,:,i], index = rowData.metab_id, columns = rowData.metab_id, dtype = float) for i in range(prec_array.shape[2])}

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





# def getComponents(adj):
#     def _dfs(adj, visited, s, res):
#         row = adj.iloc[s]
#         name = row.name
#         visited.append(name)
#         res.append(name)

#         for i in range(row.shape[0]):
#             if row.iloc[i] and adj.index[i] not in visited:
#                 _dfs(adj, visited, i, res)
#     visited = []
#     components = []
#     for s in range(adj.shape[0]):
#         metab = adj.index[s]
#         if metab not in visited:
#             res = []
#             _dfs(adj, visited, s, res)
#             components.append(res)
#     return components

# components = getComponents(adj_form)
# components.__len__()

if 1: 
    print("True")
# \\\\
# \\\\
# Save COmponent metab and pathways
# \\\\
# \\\\

def format_metab(metab_id):
    # Get pathway info for these metab_ids
    mask = rowData.metab_id.eq(metab_id)
    chem = rowData.chemical_name[mask].values[0]
    super_pw = rowData.super_pathway[mask].values[0]
    sub_pw = rowData.sub_pathway[mask].values[0]
    sep = '•'
    metab_str = f"{super_pw}{sep}{sub_pw}{sep}{chem}"
    return metab_str

def save_metab_pathways(metab_ids, out_dir, out_nm, mode='w'):
    """Save metab pathways to text file"""
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    metab_ids = list(metab_ids)
    # Get pathway info for these metab_ids
    lines = [format_metab(x) for x in metab_ids]
    line = "||".join(lines)
    # Write to file
    with open(out_dir / out_nm, mode, encoding="UTF8") as f:
        print(line, file=f)

def save_top_comps(G, out_dir, out_nm, N):
    """Find top N components and save metab pathways for each on separate lines"""
    # Get connected components sorted by size
    conn_comp = list(nx.connected_components(G))
    lngs = [len(x) for x in conn_comp]
    top_N_indices = np.argsort(lngs)[-N:]
    # Save first component (overwrites file)
    save_metab_pathways(list(conn_comp[top_N_indices[0]]), out_dir, out_nm, mode='w')
    # Append remaining components
    for idx in top_N_indices[1:]:
        save_metab_pathways(list(conn_comp[idx]), out_dir, out_nm, mode='a')

save_top_comps(G = G_form,out_dir = 'results/005/form-metabs/', out_nm = 'comp-metabs.txt',N = 10)
save_top_comps(G = G_curr, out_dir = 'results/005/curr-metabs/', out_nm = "comp-metabs.txt" , N = 10)



# \\\\
# \\\\
# Look at node overlaps
# \\\\
# \\\\


# metabolites kept after truncation/filtering
inds_curr, inds_form = [x.columns.values for x in (adj_curr, adj_form)]

# nodes unique to current
inds_distinct_curr = np.setdiff1d(inds_curr, inds_form)
# nodes unique to former
inds_distinct_form= np.setdiff1d(inds_form, inds_curr)

save_metab_pathways(metab_ids=inds_distinct_curr,out_dir= 'results/005/curr-metabs',out_nm= "unique-nodes.txt")
save_metab_pathways(metab_ids=inds_distinct_form, out_dir='results/005/form-metabs', out_nm="unique-nodes.txt")

# all nodes
inds_union = np.union1d(*(inds_curr, inds_form))
# where current has node in all nodes
in_union_curr = np.isin(inds_union, inds_curr)
# where former has node in all nodes
in_union_form  = np.isin(inds_union, inds_form)

inds_distinct_curr

def jaccard(s1, s2):
    """Intersection/union"""
    return  len(np.intersect1d(s1,s2))/len(np.union1d(s1,s2))
     
from upsetplot import UpSet, from_indicators
upset_data = from_indicators(["former", "current"],
                             pd.DataFrame({"former": in_union_form,"current": in_union_curr,},
                                          index=inds_union)
                            )
up = UpSet(upset_data.copy(), subset_size="count", element_size=60)
fig = plt.figure(figsize = (20,20))
up.plot(fig = fig)
plt.title(f"Node Jaccard= {jaccard(*(inds_curr, inds_form)).__round__(2)}\n for tau = {tau}\nand/or top N = {top_N}")
plt.savefig('results/005/node-upset.pdf')



# \\\\
# \\\\
# Look at edge overlaps for shared nodes
# \\\\
# \\\\

# shared metabolites
inds_intersect = np.intersect1d(*(inds_curr, inds_form))

# adjacencies at shared nodes
adj_shared_curr = adj_curr.loc[inds_intersect, inds_intersect].astype(bool)
adj_shared_form = adj_form.loc[inds_intersect, inds_intersect].astype(bool)

adj_shared_1d_curr = adj_shared_curr.values.flatten()
adj_shared_1d_form = adj_shared_form.values.flatten()
adj_shared_1d_union = adj_shared_1d_curr | adj_shared_1d_form

edge_jaccard = sum(adj_shared_1d_curr & adj_shared_1d_form) / sum(adj_shared_1d_curr | adj_shared_1d_form)

upset_data = from_indicators(['current', 'former'], pd.DataFrame({
    'current':adj_shared_1d_curr[adj_shared_1d_union], 
    'former':adj_shared_1d_form[adj_shared_1d_union],
}))

up = UpSet(upset_data.copy(), subset_size="count", element_size=60)
fig = plt.figure(figsize=(20, 20))
up.plot(fig=fig)
plt.title(f"Edge Jaccard = {round(edge_jaccard, 3)}\nfor tau = {tau}\nand/or top N = {top_N}")
plt.savefig('results/005/edge-upset.pdf')


#\\\\\
#\\\\\
# Distinct Edges for shared nodes
#\\\\\
#\\\\\
# edges unique to group where nodes shared
adj_edge_distinct_curr = adj_shared_curr &  ~adj_shared_form
adj_edge_distinct_form = ~adj_shared_curr &  adj_shared_form

def get_edge_node_pairs(adj):
    """get node pairs for edges"""
    return [(inds_intersect[x], inds_intersect[y]) for x,y in  zip(*np.where(adj)) ]

edge_distinct_curr = get_edge_node_pairs(adj_edge_distinct_curr)
edge_distinct_form = get_edge_node_pairs(adj_edge_distinct_form)

print(f"{int(len(edge_distinct_form)/2)} unique edges in former smokers\n{int(len(edge_distinct_curr)/2)} unique edges in current smokers")

# chem_edge_distinct_curr = [(rowData.chemical_name.values[rowData.metab_id == x][0], rowData.chemical_name.values[rowData.metab_id == y][0]) for x,y in edge_distinct_curr]
# chem_edge_distinct_form = [(rowData.chemical_name.values[rowData.metab_id == x][0], rowData.chemical_name.values[rowData.metab_id == y][0]) for x,y in edge_distinct_form]

def filter_edges(edge_list):
    """filter for where not repeat but rev order"""
    filt = []
    seen = set()
    for x,y in edge_list:
        edge = tuple(sorted([x,y]))
        if edge not in seen:
            seen.add(edge)
            filt.append(edge)
    return filt

with open( "results/005/curr-metabs/unique-edges-at-shared-nodes.txt", 'wt') as f:
    for x, y in filter_edges(edge_distinct_curr):
        x,y = (format_metab(x),format_metab(y))
        f.write(f"{x}\t{y}\n")

with open( "results/005/form-metabs/unique-edges-at-shared-nodes.txt", 'wt') as f:
    for x, y in filter_edges(edge_distinct_form):
        x,y = (format_metab(x),format_metab(y))
        f.write(f"{x}\t{y}\n")
