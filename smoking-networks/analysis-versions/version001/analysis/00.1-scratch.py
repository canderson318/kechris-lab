#!/Users/canderson/miniconda3/envs/smoknet-env/bin/python

import numpy as np
from sklearn.preprocessing import Normalizer 

x = np.linspace(0, 20, num=100)

print(x.shape)

x_n = Normalizer('l2').fit_transform(x.reshape(-1,1))

print(x_n.shape)
print(x_n)


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
