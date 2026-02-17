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
##
# filter warning
warnings.filterwarnings('ignore', message='.*chained assignment.*')
##
# assign highest level directory
op_sys = platform.system()

# set root depending on system
if op_sys == 'Linux':
    root=Path.home() / '/projects/canderson2@xsede.org/kechris-lab/smoking-networks/'
elif op_sys == "Darwin": 
    root=Path.home() / 'Documents/school/local-kechris-lab/kechris-lab/smoking-networks'
else:
    raise SystemError("System Not Identifiable")

# try changint to root
try:
    os.chdir(root / 'analysis-versions/version001')
except OSError as e:
    print(f"Failed to set working directory: {e}", file=sys.stderr)
    sys.exit(1)
##
# set results directory
Path("results/004").mkdir(exist_ok=True)

# \\\\
# \\\\
# Load Data 
# \\\\
# \\\\

rowDat = pd.read_csv('processed-data/001/rowData.csv', index_col = 0)
rcfgl = pkl.load(open('results/003/RCFGL-output/RCFGL.pkl','rb'))

len(rcfgl)
prec_array, _, _ = rcfgl
prec_array.shape

# make precisions dataframes
nms = ['curr', "form"]
prec_dict = {nms[i]: pd.DataFrame(prec_array[:,:,i], index = rowDat.metab_id, columns = rowDat.metab_id, dtype = float) for i in range(prec_array.shape[2])}


# \\\\
# \\\\
# Process Network
# \\\\
# \\\\

#\\
## make adjacency
#\\

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
theta = list(prec_dict.values())[0]
truncation_value = tau = 0.04
top_N = "all"

G_curr, adj_curr,  G_form, adj_form = adjGraphs(list(prec_dict.values()), tau, top_N)

# check my function is correct
adj_formx, _ = MakeAdjMatrix(theta = prec_dict['form'], truncation_value=tau, top_N = top_N, names = 'none')
adj_formx.shape
adj_form.shape
if not np.array_equal( adj_form.to_numpy(), adj_formx) :
    raise ValueError("My adjacency function returns a value inequivalent to original function.")

# \\\\
# \\\\
# Calculate network metrics
# \\\\
# \\\\


# \\
# look at how different tau changes num edges and their strength
# \\

taus = np.arange(0.001, 0.05, step = 0.005)
result = {'curr':[], 'form':[]}  # Start with empty lists
for tau in taus:
    G_curr, adj_curr, G_form, adj_form = adjGraphs(list(prec_dict.values()), tau, top_N)
    n_edge_curr = (adj_curr > 0).sum().sum()
    n_edge_form = (adj_form > 0).sum().sum()
    
    # Append values to dictionary lists
    result['curr'].append(n_edge_curr)
    result['form'].append(n_edge_form)

result = pd.DataFrame(result)
result['tau'] = taus
result = result.melt(id_vars="tau")

# sns.lineplot(result, x = 'tau', y = 'value', hue = 'variable')
# plt.show()


#\\
# look at edge weight distributions
#\\

# pseudo partial correlations
partial_curr= CovtoCor(prec_dict['curr'])
partial_form = CovtoCor(prec_dict['form'])

form, curr = (partial_form.values.flatten() ,  partial_curr.values.flatten())

fig, ax = plt.subplots(figsize=(8, 4))
sns.kdeplot(form, ax=ax, color="steelblue", label="form",alpha = .5,clip = (None, 1) )
sns.kdeplot(curr, ax=ax, color="orange", label="curr", alpha = .5, clip = (None, 1) )
ax.legend()
plt.title(f"Distributions of partial-correlation (tau = {tau})")
plt.savefig('results/004/partial-corr-distr.pdf')
plt.close()

# transform to less skewed for plotting
log_form, log_curr = [np.log(abs(x)) for x in [form,curr]]

fig, ax = plt.subplots(figsize=(8, 4))
sns.kdeplot(log_form, ax=ax, color="steelblue", label="form",alpha = .9,clip = (None, 0) )
sns.kdeplot(log_curr, ax=ax, color="orange", label="curr", alpha = .9, clip = (None, 0) )
ax.legend()
plt.title(f"Distributions of log(abs(partial-correlation)>0) (tau = {tau})")
plt.savefig('results/004/log-partial-corr-distr.pdf')
plt.close()


df = pd.DataFrame({"wgt_form": form,  "wgt_curr": curr})
fig, ax = plt.subplots(figsize = (8,8))
sns.scatterplot(data =df, x = "wgt_curr", y = "wgt_form")
plt.xlabel("Current Smoker Weights")
plt.ylabel("Former Smoker Weights")
ax.set_title("Former versus Current Smoker Edge Weights (Pseudo Partial Correlation)")
plt.savefig("results/004/former-current-weight-scatter.png")
plt.close()

log_df = df.apply(np.abs).apply(np.log).dropna()

fig, ax = plt.subplots(figsize = (8,8))
sns.scatterplot(data =log_df, x = "wgt_curr", y = "wgt_form")
plt.xlabel("Current Smoker Weights")
plt.ylabel("Former Smoker Weights")
ax.set_title("Former versus Current Smoker Edge Weights (log[abs(partial-correlation) > 0])")
plt.savefig("results/004/log-former-current-weight-scatter.png")
plt.close()

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
df.index = df.iloc[:,0]

# add sum edges >/< 0 
new_row = pd.DataFrame({df.columns[0]: ["Sum Negative Edges"], 'Former': [(prec_dict['form'] < 0).sum().sum()], 'Current': [(prec_dict['curr'] < 0).sum().sum()] }, index=["Sum Negative Edges"])
df = pd.concat([df, new_row])
new_row = pd.DataFrame({df.columns[0]: ["Sum Positive Edges"], 'Former': [(prec_dict['form'] > 0).sum().sum()], 'Current': [(prec_dict['curr'] > 0).sum().sum()] }, index=["Sum Positive Edges"])
df = pd.concat([df, new_row])

# save to csv
df.to_csv(f'results/004/tau{tau}-network-metrics.csv', index = False)

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

df = pd.DataFrame({"degree": np.concatenate([deg_form, deg_curr]),"group": (["form"] * len(deg_form)) + (["curr"] * len(deg_curr)),})

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

df = pd.DataFrame({"length": np.concatenate([short_paths_form, short_paths_curr]), 'group': (["form"] * len(short_paths_form)) + (["curr"] * len(short_paths_curr))})

fig, ax = plt.subplots(figsize=(8, 5))
sns.histplot(data=df,x="length",hue="group",fill=True,alpha=0.5,common_norm=False,multiple="dodge",  ax=ax, kde=True,kde_kws={"bw_adjust":2})
# sns.kdeplot(data=df,x="length",hue="group",fill=True,alpha=0.5,common_norm=False,  ax=ax)
ax.set_title("Shortest Path Lengths")
ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
plt.savefig("results/004/shortest-path-densities.pdf")


