### MAKE DATA DEMOGRAPHIC TABLE ###
### MAKE DATA DEMOGRAPHIC TABLE ###
### MAKE DATA DEMOGRAPHIC TABLE ###
### MAKE DATA DEMOGRAPHIC TABLE ###

import pandas as pd
import numpy as np
import sys
from pathlib import Path
import os
import platform
import subprocess as sp



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

# make results/processed-data dirs
Path("results/002_3").mkdir(exist_ok=True)


# \\\\
# \\\\
# Load Data
# \\\\
# \\\\

rowData = pd.read_csv('processed-data/001/rowData.csv', index_col=0)
colData = pd.read_csv('processed-data/001/colData.csv', index_col=0)
counts_curr = np.loadtxt('processed-data/002/separate/current.csv', delimiter=',')
counts_form = np.loadtxt('processed-data/002/separate/former.csv', delimiter=',')

counts_curr.shape
counts_form.shape


#\\\
#\\\
# Format Data
#\\\
#\\\

# separate into smoking/non-smoking matrices
# smoking_status	0 = Never Smoker; 1 = Former smoker; 2 = Current smoker			
nev_smok_mask = colData.smoking_status == 0
form_smok_mask = colData.smoking_status == 1
curr_smok_mask = colData.smoking_status == 2

print(f"\n{nev_smok_mask.sum()} never smoked\n{form_smok_mask.sum()} smoked in the past\n{curr_smok_mask.sum()} currently smoke\n")

# subset for people i am using
sub_colData = colData.loc[form_smok_mask | curr_smok_mask, :]


# \\\
# \\\
# Tables
# \\\
# \\\


cols = ['age_visit', 'bmi']
vars = ["age", "BMI"]
sexf = sub_colData.sexf
smok = sub_colData.smoking_status

with open("results/002_3/dem-tab.csv", 'wt') as f:
    print(f"Var,SexF,Mean,Median,Std,Current Smokers,Former Smokers,Total", file = f)
    for sex in [True, False]:
        for v,c in zip(vars,cols):
            x = sub_colData.loc[sexf == sex, c]
            mean = x.mean().round(3)
            median = x.median().round(3)
            std = x.std().round(3)
            if c == "bmi":
                counts = smok.loc[sexf == sex].value_counts().values
                print(f"{v},{sex},{mean},{median},{std},{counts[0]},{counts[1]},{counts[0]+counts[1]}", file = f)     
            else:
                print(f"{v},{sex},{mean},{median},{std},,,", file = f)     
                
    print(f",,,,,{counts_form.shape[0]},{counts_curr.shape[0]},{sub_colData.shape[0]}", file = f)

# metab info
with open("results/002_3/metab-info.txt", 'wt') as f:
    print(f"{rowData.shape[0]}", file = f)
    print(", ".join(rowData.super_pathway.unique()), file = f) 

# \\\\
# \\\\
# Save them
# \\\\
# \\\\
def cp(path):
    sp.run(["cp", path, 
            '/Users/canderson/Documents/school/local-kechris-lab/presentations/rotation-pres/data'])

cp('results/002_3/dem-tab.csv')
cp('results/002_3/metab-info.txt')
