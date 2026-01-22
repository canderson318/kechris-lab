#!/Users/canderson/miniconda3/envs/smoknet-env/bin/python

import pandas as pd
import numpy as np
import sys
import os
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import seaborn as sns
import plotly.graph_objects as go
import kaleido

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
        # method='linear'   # 'nearest' 'linear' 'cubic'
        method='cubic'
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


# interactive plot
df, col_label, title, filename = (
    res.iloc[(res.l1.values<.06) & (res.l2.values <.04) , :], 
    '√AIC',
    'Interpolated AIC surface (zoomed)',
    'results/002/interactive-2x-zoomed-AIC-grid.html'
)

# interpolation grid
l1_grid = np.linspace(df.l1.min(), df.l1.max(), 500)
l2_grid = np.linspace(df.l2.min(), df.l2.max(), 500)
L1, L2 = np.meshgrid(l1_grid, l2_grid)

AIC_grid = griddata(
    points=(df.l1, df.l2),
    values=np.sqrt(df.AIC),
    xi=(L1, L2),
    # method='linear'
    method='cubic'
)

fig = go.Figure()

# interpolated surface
fig.add_trace(
    go.Heatmap(
        x=l1_grid,
        y=l2_grid,
        z=AIC_grid,
        colorscale="Viridis",
        colorbar=dict(title=col_label),
        zsmooth="best"
    )
)

# original sample points
fig.add_trace(
    go.Scatter(
        x=df.l1,
        y=df.l2,
        mode="markers",
        marker=dict(
            size=6,
            color=np.sqrt(df.AIC),
            colorscale="Viridis",
            showscale=False,
            line=dict(width=0.5, color="black")
        ),
        name="Observed points",
        hovertemplate="l1=%{x}<br>l2=%{y}<br>√AIC=%{marker.color:.3f}<extra></extra>"
    )
)

fig.update_layout(
    title=title,
    xaxis_title="l1",
    yaxis_title="l2",
    template="plotly_white"
)

if filename is not None:
    fig.write_html(filename)


