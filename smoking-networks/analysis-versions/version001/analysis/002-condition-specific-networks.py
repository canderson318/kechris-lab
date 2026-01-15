#!/Users/canderson/miniconda3/envs/smoknet-env/bin/python

import pandas as pd
import sys
import os
from sklearn.preprocessing import StandardScaler

RCFGL_path = '/Users/canderson/Documents/school/local-kechris-lab/kechris-lab/smoking-networks/RCFGL'
os.chdir(RCFGL_path)
sys.path.insert(0, 'Python_functions')
from RCFGL import RFGL, RCFGL

os.chdir(RCFGL_path)
from Dstream_functions import*


