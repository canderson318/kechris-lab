
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
import re
from itertools import chain

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
out_dir = Path("results/005_3")
out_dir.mkdir(exist_ok=True)

#\\\\
#\\\\
# Load Data
#\\\\
#\\\\
rowData = pd.read_csv('processed-data/001/rowData.csv', index_col = 0)
rowData.set_index("metab_id", inplace=True)
colData = pd.read_csv('processed-data/001/colData.csv', index_col=0)
adj_curr, adj_form = [pd.read_csv(p, index_col= 0) for p in ['results/005/curr-adj.csv','results/005/form-adj.csv']]

from analysis.utils.myDstream_functions import MetabLoader

comp_curr = MetabLoader("results/005/curr-metabs/comp-metabs.txt")
metab_comp_curr = comp_curr.load(at = 2, rowData=rowData, id = False)
comp_form = MetabLoader("results/005/form-metabs/comp-metabs.txt")
metab_comp_form = comp_form.load(at = 2, rowData=rowData, id = False)

#\\\\
#\\\\
# Compare Sphingomyelin sizes
#\\\\
#\\\\

# component metabolites
SM_curr = metab_comp_curr[5]
SM_form = metab_comp_form[5]

# # unique metabolites
# SM_curr = list(rowData.loc[np.setdiff1d(adj_curr.index, adj_form.index), 'chemical_name'])
# SM_form = list(rowData.loc[np.setdiff1d(adj_form.index, adj_curr.index), 'chemical_name'])

# all metabolites
all_SM_curr = list(rowData.loc[adj_curr.index, 'chemical_name'])
all_SM_form = list(rowData.loc[adj_form.index, 'chemical_name'])

# filter for sphingolipids
SM_curr = [x for x in SM_curr if "sphingomyelin" in x]
SM_form = [x for x in SM_form if "sphingomyelin" in x]
all_SM_curr = [x for x in all_SM_curr if "sphingomyelin" in x]
all_SM_form = [x for x in all_SM_form if "sphingomyelin" in x]



# take avg of each sm d#:#/#:# for both isobarics
# compare distr between groups
# compare number connected between each type/size of sm
# plot FA size against degree

def extract_nums(text):
    """Extract number groups from sphingomyelin notation. Example: (d18:1/20:0, d16:1/22:0) -> [(18,1,20,0), (16,1,22,0)]"""
    # Pattern to match d##:##/##:## format
    pattern = r'd(\d+):(\d+)/(\d+):(\d+)'
    # Find all matches
    matches = re.findall(pattern, text)
    # Convert strings to integers and return as tuples
    return [tuple(map(int, match)) for match in matches]


# # average isobaric species at each location
# process = lambda tuples_list: [np.array(x).mean(axis = 0) for x in tuples_list]
# nums_SM_curr = np.array(process([extract_nums(x) for x in SM_curr])).round(3)
# nums_SM_form = np.array(process([extract_nums(x) for x in SM_form])).round(3)

#  make each isobaric species own item in list
nums_SM_curr = np.array(list(chain.from_iterable([extract_nums(x) for x in SM_curr]))).round(3)
nums_SM_form = np.array(list(chain.from_iterable([extract_nums(x) for x in SM_form]))).round(3)
nums_all_SM_curr = np.array(list(chain.from_iterable([extract_nums(x) for x in all_SM_curr]))).round(3)
nums_all_SM_form = np.array(list(chain.from_iterable([extract_nums(x) for x in all_SM_form]))).round(3)

# #  choose first of isobaric species 
# nums_SM_curr = np.array([extract_nums(x)[0] for x in SM_curr]).round(3)
# nums_SM_form = np.array([extract_nums(x)[0] for x in SM_form]).round(3)
# nums_all_SM_curr = np.array([extract_nums(x)[0] for x in all_SM_curr]).round(3)
# nums_all_SM_form = np.array([extract_nums(x)[0] for x in all_SM_form]).round(3)

#\\\\
#\\\\
# plot
#\\\\
#\\\\

labels = {
    (0, 0): ("Base: #Carbons", 0),
    (1, 0): ("FA-chain: #Carbons", 2),
    (0, 1): ("Base: #Double bonds", 1),
    (1, 1): ("FA-chain: #Double bonds", 3)
}

# plot each loc distribution for each group; SM components
fig, ax = plt.subplots(2,2, figsize = (10,10))
for i in range(2):
    for j in range(2):
        label, idx = labels[(i, j)]
        
        df = pd.DataFrame({
            "group": (["curr"] * nums_SM_curr.shape[0]) + (["form"] * nums_SM_form.shape[0]),
            "x": np.concatenate([nums_SM_curr[:, idx], nums_SM_form[:, idx]])
        })
        
        sns.kdeplot(df, x="x", hue="group", fill=True, alpha=.2, ax=ax[i, j], 
                   clip=(0, None), common_norm=False)
        ax[i, j].set_xlabel(label)
        ax[i, j].set_ylabel("Density")
fig.suptitle(f"Component Distribution of Sphingomyelin Base/Chain Lengths\nNcurr = {nums_SM_curr.shape[0]}, Nform = {nums_SM_form.shape[0]}\n(all isobaric species)", 
             fontsize=16, fontweight="bold")
plt.tight_layout()
plt.savefig("results/005_3/component-sphingomyelin-base-chain-length-distrs.pdf")


# all SMS
fig, ax = plt.subplots(2,2, figsize = (10,10))
for i in range(2):
    for j in range(2):
        label, idx = labels[(i, j)]
        
        df = pd.DataFrame({
            "group": (["curr"] * nums_all_SM_curr.shape[0]) + (["form"] * nums_all_SM_form.shape[0]),
            "x": np.concatenate([nums_all_SM_curr[:, idx], nums_all_SM_form[:, idx]])
        })
        
        sns.kdeplot(df, x="x", hue="group", fill=True, alpha=.2, ax=ax[i, j], 
                   clip=(0, None), common_norm=False)
        ax[i, j].set_xlabel(label)
        ax[i, j].set_ylabel("Density")
fig.suptitle(f"Distribution of All Sphingomyelin Base/Chain Lengths\nNcurr = {nums_all_SM_curr.shape[0]}, Nform = {nums_all_SM_form.shape[0]}\n(all isobaric species)", 
             fontsize=16, fontweight="bold")
plt.tight_layout()
plt.savefig("results/005_3/sphingomyelin-base-chain-length-distrs.pdf")

#\\\
# Summaries
#\\\

def summSM(arr):
    arr = arr.round(3)
    print(f"{arr[0]} base carbons, {arr[1]} base double-bonds,\n{arr[2]} FA-chain carbons, and {arr[3]} FA-chain double-bonds")

from scipy import stats

f = lambda arr: arr.mean(axis = 0)
print("Former:")
summSM(f(nums_SM_form))
print("Current:")
summSM(f(nums_SM_curr))

rowData.columns