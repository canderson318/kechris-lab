# Copied from ../../../../RCFGL/Python_functions/Dstream_functions.py

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import networkx as nx

# My versions of adjacnency functions 
#\\\\#\\\\#\\\\#\\\\#\\\\#\\\\#\\\\#\\\\#\\\\
def Adjacency(theta, truncation_value, top_N):
    """Make an adjacency matrix for where edges > truncation value"""
    ## save names
    # names = np.arange(theta.shape[0])
    indx = theta.index
    # calculate partial correlation strength
    ppcorr = np.abs(CovtoCor(theta))
    # fill with zeros where below truncation
    ppcorr[ppcorr < truncation_value] = 0
    # fill lower triangle including diagonal with zeros (undirected graph so upper and lower are the same)
    uppr = np.triu(ppcorr, 1)
    # number upper tri cells where non-zero
    if top_N == "all":
        top_N = len(np.where(uppr != 0)[0])
    # find the indices of the top N largest values in your flattened array:
    idx = np.argpartition(uppr, uppr.size - top_N, axis=None)[-top_N:]
    # Convert 1D cell locations to 2d array indices
    locations = np.column_stack(np.unravel_index(idx, uppr.shape))
    # set cells at locations to 1
    adj = np.zeros(uppr.shape)
    for k in locations:
        adj[k[0], k[1]] = 1
        adj[k[1], k[0]] = 1
    #
    # get nonzero columns
    nonz = np.where(np.sum(adj, axis = 0) != 0)[0]
    # nonzero names
    nonz_names = indx[nonz]
    # remove columns/rows with only zero entries
    adj = adj[np.ix_(nonz, nonz)]
    return pd.DataFrame(adj, index=nonz_names, columns=nonz_names, dtype = int)

def CovtoCor(covariance):
    """Calculate (pseudo) Partial Correlation: converts precision to dependent correlation strengths"""
    v = np.sqrt(np.diag(covariance))
    outer_v = np.outer(v, v)
    correlation = covariance / outer_v
    return correlation



class MetabLoader:
    """Loader for metabolite data from formatted text files."""
    
    def __init__(self, path):
        self.path = path
        self._lines = None
    
    def _read_lines(self):
        """Read and cache lines from the file."""
        if self._lines is None:
            with open(self.path, 'r') as f:
                self._lines = f.readlines()
        return self._lines
    
    def _parse_components(self, at):
        lines = self._read_lines()
        return [[y.split("•")[at] for y in x.split("||")] for x in lines]
    
    def _get_metab_ids(self, comp_list, rowData):
        return [
            rowData.metab_id.values[rowData.chemical_name.isin(component)] 
            for component in comp_list
        ]
    
    def load(self, at, id=True, rowData=None):
        metabs = self._parse_components(at)
        
        if id:
            if rowData is None:
                raise ValueError("rowData required when id=True")
            return self._get_metab_ids(metabs, rowData)
        
        return metabs


#\\\\#\\\\#\\\\#\\\\#\\\\#\\\\#\\\\#\\\\#\\\\


# OLD
def MakeAdjMatrix(theta, truncation_value, top_N, names):
    theta = abs(CovtoCor(theta)); 
    if type(names) == str:
       names = np.array(range(0, theta.shape[0]))
    else:
        names = np.array(names)
    arr = np.triu(theta, 1)
    theta[theta < truncation_value] = 0
    if top_N == 'all':
       theta[theta < truncation_value] = 0 
       arr = np.triu(theta, 1)
       top_N = len(np.where(arr != 0)[0])

    idx = np.argpartition(arr, arr.size - top_N, axis=None)[-top_N:]
    locations = np.column_stack(np.unravel_index(idx, arr.shape))
    adj = np.zeros(arr.shape)
    for k in locations:
     adj[k[0], k[1]] = adj[k[1], k[0]] = 1
    
    nonz = np.where(np.sum(adj, axis = 0) != 0)
    nonz_names = names[nonz[0]]
    adj = adj[np.ix_(nonz[0],nonz[0])]
    Adjacency = []
    Adjacency.append(adj)
    Adjacency.append(nonz_names)
    return Adjacency

