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

# root=Path.home() / 'Documents/school/local-kechris-lab/kechris-lab/smoking-networks'
root= Path('/projects/canderson2@xsede.org/kechris-lab/smoking-networks')
RCFGL_path = root / 'RCFGL'
os.chdir(RCFGL_path)
try: 
    sys.path.remove("Python_functions")
except:
    print('Python_functions not in $PATH')

sys.path.insert(0, 'Python_functions')
from RCFGL import RFGL, RCFGL
from Dstream_functions import*

os.chdir(root / 'analysis-versions/version001')


#\\\
#\\\
#  Load Data
#\\\
#\\\

nms=  ['current', 'former']

all = []
for nm in nms :
    all.append( np.loadtxt(f"processed-data/002/separate/{nm}.csv", delimiter = ',', ) )

if len(all)<2:
    raise ValueError("Precision matrices not loaded")

[x.shape for x in all]

#\\\
#\\\
# Run RCFGL
#\\\
#\\\
l1 = 0.316
l2 = 0.158

print(f"Running RCFGL with lambda1 = {l1} and lambda2 = {l2}")
RCFGL_output = RCFGL(A = all, 
                   ADMMmaxiter = 100, 
                   admmtol = 0.001,
                   lambda1 = l1, 
                   lambda2 = l2) # https://doi.org/10.1371/journal.pcbi.1010758
print("RCFGL Finished")

precision_matrices_array, AIC, time = RCFGL_output

print(f"AIC = {round(AIC[0],3)}, runtime = {time}") 
precision_matrices_array.shape

Adjacency_all = MakeAdjMatrix_all(RCFGL_output, truncation_value = 0.05, top_N = 75, names = 'default')

#\\\
#\\\
#  Save Results
#\\\
#\\\

print("Saving RCFGL Outputs...")

# make paths
precision_path = Path('results/002/RCFGL-output/precision')
precision_path.mkdir(parents=True, exist_ok=True)

adjacency_path = Path('results/002/RCFGL-output/adjacency')
adjacency_path.mkdir(parents=True, exist_ok=True)


# save arrays
for i in range(0,len(Adjacency_all)):
    np.savetxt(precision_path / f'{nms[i]}.csv', precision_matrices_array[:, :, i], delimiter = ',')
    np.savetxt(adjacency_path / f'{nms[i]}.csv', Adjacency_all[i], delimiter = ',')

print("RCFGL Saved")

