#!/Users/canderson/miniconda3/envs/smoknet-env/bin/python

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import seaborn as sns

# assign highest level directory
root=Path.home() / 'Documents/school/local-kechris-lab/kechris-lab/smoking-networks'
os.chdir(root / 'analysis-versions/version001')

#\\\
#\\\
# Load in AICs
#\\\
#\\\

files = os.listdir('results/002/aics')

res = pd.DataFrame([np.loadtxt(f'results/002/aics/{f}', delimiter = ',') for f in files],  columns = ['l1', 'l2', 'AIC'])


#\\\
#\\\
# Plot AIC across l1/l2
#\\\
#\\\

plt.scatter(res.l1, res.l2, c = res.AIC)

l1_grid = np.linspace(res.l1.min(), res.l1.max(), 300)
l2_grid = np.linspace(res.l2.min(), res.l2.max(), 300)

L1, L2 = np.meshgrid(l1_grid, l2_grid)

AIC_grid = griddata(
    points=(res.l1, res.l2),
    values=np.sqrt(res.AIC),
    xi=(L1, L2),
    method='linear'   # 'nearest' or 'cubic' also possible
)


plt.figure(figsize=(10, 9))
plt.imshow(
    AIC_grid,
    origin='lower',
    extent=(l1_grid.min(), l1_grid.max(), l2_grid.min(), l2_grid.max()),
    aspect='auto'
)
plt.colorbar(label='AIC')

plt.xlabel('l1')
plt.ylabel('l2')
plt.title('Interpolated AIC surface')
plt.show()

