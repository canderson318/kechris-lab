#!/Users/canderson/miniconda3/envs/smoknet-env/bin/python

import pandas as pd
import numpy as np
import sys
import os
from numpy.typing import NDArray
from typing import Any
import re
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.preprocessing import StandardScaler

root=Path.home() / 'Documents/school/local-kechris-lab/kechris-lab/smoking-networks'

os.chdir(root / 'analysis-versions/version001')

# make results/processed-data dirs
for p in ["results/002", 'processed-data/002'] :
    Path(p).mkdir(exist_ok=True)


#\\\
#\\\
#  Load Data
#\\\
#\\\

# os.listdir('processed-data/001')
rowData = pd.read_csv('processed-data/001/rowData.csv', index_col=0)
colData = pd.read_csv('processed-data/001/colData.csv', index_col=0)
# counts = pd.read_csv('processed-data/001/raw_counts.csv', index_col=0)
adjusted_logcounts = pd.read_csv('processed-data/001/adjusted_logcounts.csv', index_col=0)

# transpose to subj x metabolites
adjusted_logcounts = adjusted_logcounts.T

#\\\
#\\\
#  Format Data
#\\\
#\\\

# separate into smoking/non-smoking matrices
# smoking_status	0 = Never Smoker; 1 = Former smoker; 2 = Current smoker			
nev_smok_mask = colData.smoking_status == 0
form_smok_mask = colData.smoking_status == 1
curr_smok_mask = colData.smoking_status == 2

print(f"{nev_smok_mask.sum()} never smoked\n{form_smok_mask.sum()} smoked in the past\n{curr_smok_mask.sum()} currently smoke")

nev_smok_logcounts = adjusted_logcounts.loc[nev_smok_mask , :]
curr_smok_logcounts = adjusted_logcounts.loc[curr_smok_mask , :]
form_smok_logcounts = adjusted_logcounts.loc[form_smok_mask , :]

# Scaling the dataframes such that the columns have 0 mean
def process(df: pd.DataFrame)->pd.DataFrame:
    scaler = StandardScaler(with_std=False)
    df_scaled = scaler.fit_transform(df)
    df_scaled = pd.DataFrame(df_scaled, index = df.index, columns = df.columns)
    return df_scaled

nev_scaled = process(nev_smok_logcounts)
form_scaled = process(form_smok_logcounts)
curr_scaled = process(curr_smok_logcounts)

# Combined list of the scaled dataframes
all = [np.array(nev_scaled), np.array(form_scaled), np.array(curr_scaled)]

nms = ['never', 'former','current']
all_dict = {nms[index]: value for index, value in enumerate(all)}

# make directory for scaled counts matrices separated by smoking status
out_path = Path('processed-data/002/separate-scaled')
out_path.mkdir(exist_ok=True)

for nm in nms:
    print(f"Saving `{nm}.csv` to {out_path}")
    np.savetxt(Path('processed-data/002/separate-scaled') / f'{nm}.csv',  all_dict[nm], delimiter = ',')
