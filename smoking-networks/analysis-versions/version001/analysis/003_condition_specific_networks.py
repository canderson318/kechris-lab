### CONDITION SPECIFIC NETWORKS ###
### CONDITION SPECIFIC NETWORKS ###
### CONDITION SPECIFIC NETWORKS ###
### CONDITION SPECIFIC NETWORKS ###

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
import pickle as pkl
import platform
import warnings

warnings.filterwarnings("ignore", message='.*import pkg_resources as pr.*')

# assign highest level directory
op_sys = platform.system()

if op_sys == 'Linux':
    root=Path.home() / '/projects/canderson2@xsede.org/kechris-lab/smoking-networks/'
elif op_sys == "Darwin": 
    root=Path.home() / 'Documents/school/local-kechris-lab/kechris-lab/smoking-networks'
else:
    raise SystemError("System Not Identifiable")

RCFGL_path = root / 'RCFGL'
os.chdir(RCFGL_path)
try: 
    sys.path.remove("Python_functions")
except:
    print('Python_functions not in $PATH, adding')
    sys.path.insert(0, 'Python_functions')  

from RCFGL import  RCFGL

try:
    os.chdir(root / 'analysis-versions/version001')
except OSError as e:
    print(f"Failed to set working directory: {e}", file=sys.stderr)
    sys.exit(1)


#\\\
#\\\
#  Load Data
#\\\
#\\\

nms=  [ 'current', 'former'] # order critical for downstream

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



#\\\
#\\\
#  Save Results
#\\\
#\\\

print("Saving RCFGL...")
out_path = Path('results/002/RCFGL-output/')
out_path.mkdir(parents=True, exist_ok=True)
out_file = 'RCFGL.pkl'
with open(out_path/out_file, 'wb')  as f:
    pkl.dump(obj=RCFGL_output, file =  f)
print("RCFGL Saved")

