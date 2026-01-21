#!/Users/canderson/miniconda3/envs/smoknet-env/bin/python

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import warnings
import re

# assign highest level directory
root=Path.home() / 'Documents/school/local-kechris-lab/kechris-lab/smoking-networks'
os.chdir(root / 'analysis-versions/version001')

#\\\
#\\\
# Load in AICs
#\\\
#\\\

files = os.listdir('results/002/20-jan-aics')

l1s = [
    float(m.group(1))
    for x in files
    if (m := re.search(r"l1_?(\d+(?:\.\d+)?)", x))
]


pd.Series(l1s).value_counts()