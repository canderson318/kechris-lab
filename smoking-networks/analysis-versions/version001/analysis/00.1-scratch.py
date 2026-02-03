
import numpy as np
import sklearn.datasets 
import pandas as pd

# functions for network analysis
from analysis.utils.myDstream_functions import *

iris = sklearn.datasets.load_iris(as_frame = True)['data']

cov = pd.DataFrame(np.cov(iris.T), index = iris.columns, columns = iris.columns)
prec = pd.DataFrame(np.linalg.inv(cov), index = iris.columns, columns = iris.columns)

# marginal cor
cor = CovtoCor(cov)

# partial correlations (correlations correcting for confouding relationships)
partial_cor = np.array(np.abs(-CovtoCor(prec))); np.fill_diagonal(partial_cor, 1); partial_cor =  pd.DataFrame(partial_cor, index = iris.columns, columns = iris.columns)

adj, nms = MakeAdjMatrix(prec, .5, "all", 'default')
adj; nms

G = nx.from_numpy_array(adj, create_using=nx.MultiGraph,)
labels = dict(zip(G.nodes(),iris.columns.values[nms]))
nx.set_node_attributes(G,labels ,name="label")


fig = plt.figure(figsize=(8, 8))
pos = nx.spring_layout(G, seed=120349)
nx.draw(G, pos, labels = labels, font_size=8)
plt.savefig('results/example-net.png')
plt.close()



# import numpy as np
# from sklearn.preprocessing import Normalizer 
# x = np.linspace(0, 20, num=100)
# print(x.shape)
# x_n = Normalizer('l2').fit_transform(x.reshape(-1,1))
# print(x_n.shape)
# print(x_n)


# #//
# # plot distributions before/after
# #//
# cotin_ids = rowData.metab_id[rowData.chemical_name.str.contains("cotin")]

# positions = curr_scaled.columns.get_indexer(cotin_ids)
# positions = positions[positions != -1]

# metab_idx = positions[0]

# fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

# # ---- Current smokers ----
# x = curr_smok_logcounts.iloc[:, metab_idx]
# y = form_smok_logcounts.iloc[:, metab_idx]
 
# sns.kdeplot(x, label="Current Smokers", fill=True, clip=(None, None), ax=axes[0])
# sns.kdeplot(y, label="Former Smokers", fill=True, clip=(None, None), ax=axes[0 ])
 
# axes[0].set_title("Pre Scaling")
# axes[0].set_xlabel("Metabolite Intensity")
# axes[0].set_ylabel("Density")
# axes[0].legend()

# # ---- Former smokers ----
# x = curr_scaled.iloc[:, metab_idx] 
# y = form_scaled.iloc[:, metab_idx]

# sns.kdeplot(x, label="Current Smokers", fill=True, clip=(None, None), ax=axes[1])
# sns.kdeplot(y, label="Former Smokers", fill=True, clip=(None, None), ax=axes[1])


# axes[1].set_title("Post Scaling")
# axes[1].set_ylabel("Density")
# axes[1].set_xlabel("Metabolite Intensity")
# axes[1].legend()

# plt.tight_layout()
# plt.show()
