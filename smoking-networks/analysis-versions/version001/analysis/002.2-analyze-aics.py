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
# root=Path.home() / 'Documents/school/local-kechris-lab/kechris-lab/smoking-networks'
root=Path.home() / '/projects/canderson2@xsede.org/kechris-lab/smoking-networks/'
os.chdir(root / 'analysis-versions/version001')

#\\\
#\\\
# Load in AICs
#\\\
#\\\

aic_dir = 'results/002/22-01-2026-aics'
files = os.listdir(aic_dir)

res = pd.DataFrame([np.loadtxt(Path(aic_dir) / f, delimiter = ',') for f in files],  columns = ['l1', 'l2', 'AIC'])

# how many did it complete of grid?
grid = pd.read_csv('results/002/lambda-grid.txt', delimiter=',', header =None) 
grid.columns = ['l1', 'l2']

grid_pairs = set(map(tuple, grid[['l1', 'l2']].to_numpy()))
res_pairs  = set(map(tuple, res[['l1', 'l2']].to_numpy()))

print(f"{len(grid_pairs & res_pairs)} out of {grid.shape[0]} grid value pairs completed")


#\\\
#\\\
# Plot AIC across l1/l2
#\\\
#\\\

plt.scatter(res.l1, res.l2, c = res.AIC)

l1_grid = np.linspace(res.l1.min(), res.l1.max(), 500)
l2_grid = np.linspace(res.l2.min(), res.l2.max(), 500)

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
plt.colorbar(label='√AIC')

plt.xlabel('l1')
plt.ylabel('l2')
plt.title('Interpolated AIC surface')
plt.show()

plt.savefig('results/002/AIC-grid.png')


# look at area where AIC small
res_sub = res.iloc[(res.l1.values<.08) & (res.l2.values <.1) , :]

plt.scatter(res_sub.l1, res_sub.l2, c = res_sub.AIC)

l1_grid = np.linspace(res_sub.l1.min(), res_sub.l1.max(), 500)
l2_grid = np.linspace(res_sub.l2.min(), res_sub.l2.max(), 500)

L1, L2 = np.meshgrid(l1_grid, l2_grid)

AIC_grid = griddata(
    points=(res_sub.l1, res_sub.l2),
    values=np.sqrt(res_sub.AIC),
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
plt.colorbar(label='√AIC')

plt.xlabel('l1')
plt.ylabel('l2')
plt.title('Interpolated AIC surface (zoomed)')
plt.show()

plt.savefig('results/002/zoomed-AIC-grid.png')
