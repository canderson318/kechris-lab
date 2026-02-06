import requests
from pprint import pprint
import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import networkx as nx
import pickle as pkl
import platform
from KEGGRESTpy import kegg_get

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

# \\\\
# \\\\
# Load Data 
# \\\\
# \\\\

rowDat = pd.read_csv('processed-data/001/rowData.csv', index_col = 0)
rowDat = rowDat.iloc[:100, :]

# search compound by name
def get_top_pathway(cid):
    if pd.isna(cid):
        return None
    try:
        record = kegg_get([cid])
        if 'PATHWAY' in record.keys():
            return record['PATHWAY'][0].split('  ')[-1]
    except Exception as e:
        print(f"could not retrieve record for cid: {e}", file=sys.stderr)
    return None

pathways = [get_top_pathway(cid) for cid in rowDat.kegg.values[:10] ]
