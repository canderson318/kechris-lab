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

def plot_grid(df, col_label, title, filename):

    plt.scatter(df.l1, df.l2, c = df.AIC)

    l1_grid = np.linspace(df.l1.min(), df.l1.max(), 500)
    l2_grid = np.linspace(df.l2.min(), df.l2.max(), 500)

    L1, L2 = np.meshgrid(l1_grid, l2_grid)

    AIC_grid = griddata(
        points=(df.l1, df.l2),
        values=np.sqrt(df.AIC),
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
    plt.colorbar(label=col_label)

    plt.xlabel('l1')
    plt.ylabel('l2')
    plt.title(title)
    plt.show()

    plt.savefig(filename)

plot_grid(res,'√AIC','Interpolated AIC surface','results/002/AIC-grid.png')

# look at area where AIC small
plot_grid(res.iloc[(res.l1.values<.08) & (res.l2.values <.1) , :], '√AIC','Interpolated AIC surface (zoomed)','results/002/zoomed-AIC-grid.png')

# look at area where AIC small again
plot_grid(res.iloc[(res.l1.values<.06) & (res.l2.values <.04) , :], 
    '√AIC','Interpolated AIC surface (zoomed)','results/002/2x-zoomed-AIC-grid.png')

