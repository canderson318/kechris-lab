#!/Users/canderson/miniconda3/envs/smoknet-env/bin/python

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import warnings
from itertools import product

# suppress cppyy warning
warnings.filterwarnings(
    "ignore",
    message=r".*pkg_resources is deprecated as an API.*",
    module=r"cppyy(\..*)?$",
    category=UserWarning,
)

# assign highest level directory
# root=Path.home() / 'Documents/school/local-kechris-lab/kechris-lab/smoking-networks'
root = Path('/projects/canderson2@xsede.org/kechris-lab/smoking-networks/')

# load bespoke definitions
RCFGL_path = root / 'RCFGL'
os.chdir(RCFGL_path)#!/Users/canderson/miniconda3/envs/smoknet-env/bin/python
sys.path.insert(0, 'Python_functions')
from RCFGL import RCFGL
from Dstream_functions import*

# chang back to main project directory
os.chdir(root / 'analysis-versions/version001')


#\\\
#\\\
# Load Data
#\\\
#\\\

# load data to list
def load_dat():

    # nms=  ['never', 'current', 'former']
    nms=  ['current', 'former'] # exclude never smokers
    
    all = []
    for nm in nms :
        all.append( np.loadtxt(f"processed-data/002/separate/{nm}.csv",dtype=float, delimiter = ',') )
    
    # # |||For Testing|||
    # np.random.seed(102); inds = np.random.choice(all[0].shape[1], 100, replace = False)
    # all = [x[:,inds] for x in all ]
    # # |||^for testing|||

    return all

# # Make data shared among workers
# GLOBAL_DATA  = None

# def init_worker():
#     global GLOBAL_DATA
#     GLOBAL_DATA = load_dat()
#     print(f"[PID {os.getpid()}] Data loaded with shapes:",
#           [x.shape for x in GLOBAL_DATA])



#\\\
#\\\
# Run RCFGL
#\\\
#\\\

def run_single(params):

    l1,l2 = params

    # counts_list = GLOBAL_DATA
    counts_list = load_dat()

    _, AIC, _ = RCFGL(A = counts_list, 
                   ADMMmaxiter = 100, 
                   admmtol = 0.001,
                   lambda1 = l1, 
                   lambda2 = l2) # https://doi.org/10.1371/journal.pcbi.1010758
    
    outdir = Path("results/002/aics")
    outdir.mkdir(parents=True, exist_ok=True)

    outfile = outdir / f"AIC_l1_{l1:.3g}_l2_{l2:.3g}.txt"
    
    with open(outfile, "w") as f:
        f.write(f"{l1},{l2},{AIC[0]}\n")

    return l1, l2, AIC

# init_worker(); run_single((.1,.1)) # test

#\\\
#\\\
# Run in parallel
#\\\
#\\\

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--l1", type=float, required=True)
    parser.add_argument("--l2", type=float, required=True)
    args = parser.parse_args()

    # init_worker()
    # args = argparse.Namespace(l1=0.1, l2=0.01)

    # feed tuple of params into function as separate processes
    l1, l2, AIC = run_single((args.l1,args.l2))
    
    print(f"Finished l1={l1:.3g}, l2={l2:.3g}, AIC={AIC}")
