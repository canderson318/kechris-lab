# Smoking Metabolite Networks
Research project during my second rotation with Professor Katerina Kechris analyzing smoking effect on metabolite interactions using the [RCFGL](https://doi.org/10.1371/journal.pcbi.1010758) tool. 
Find my proposal [here](rotation-proposal.md).

#  Set up
Reference `setup.sh` for how to set up the environment. 
Key points:
- Need python/R dependencies : run `conda env create -f environment.yaml` 
- Need RCFGL library from git repo: clone Souvick's repo
    - *Note: this rcfgl version depends on older cpp compiler and other libraries so these must be installed manually outside of conda (using conda's pip though) (see `setup.sh`)*

## Bugs
The `Python_functions/get_screening_2.py` script has a deprecated `scipy.ElasticNet` parameter, `normalize`. This does not exist any more and can be removed from those lines (lines 69 & 80 I think?). Reference [analysis-versions/version001/analysis/utils/fixed_get_screening_2.py](analysis-versions/version001/analysis/utils/fixed_get_screening_2.py) for my fix. 

# Overview
I jump between python and R a lot because most of the stats/data cleaning is done using other people's scripts and RCFGL is implemented in python. 

This project has the following directory structure:

---
.
├── **analysis-versions/** (pseudo-version control) \
│   ├── **version-control.txt** (keep track of major changes that warrant a new version)\
│   └── **version001/** (only one version of analysis)\
│       ├── **analysis/** (analysis scripts named in run order )\
│       ├── **logs/** (where i store my alpine slurm logs)\
│       ├── **processed-data/** (intermediate datasets; directory number corresponds to script)\
│       ├── **results/ (outputs** including plots and summary tables; directory number corresponds to script)\
├── **info/** (notes; e.g. email that explains dataset)\
├── **raw-data/** (COPDGene metabolite counts and patient data csvs)\
├── **RCFGL/** (cloned RCFGL repo)\
├── **environment.yaml**\
└── **setup.sh**\

---

### Script info/purpose

**Extraneous Scripts**
`00_1_make_example_graphs.py`
 - make simple graph plots for illustrative purposes

`00_scratch.py`
- scratchpad script

`get_kegg.py`
- playing around with querying KEGG database

**Data Cleaning/Prep**
`001_preprocess_data.R`
- load data 
- filter out controls
- covariate adjust, log-normalize, and standardize metabolite intensities
    - uses Arun's `utils/DataPreProcessCodes.R` 

**Param Exploration**
`002_separate_conditions.py`
- make separate matrices for each condition

`002_1_runRCFGL_parallel.sh`
- bash script to run `runRCFGL.py` on different lambdas to find optimal pair with low AIC

`002_2_analyze_aics.py`
- plot AIC values on lambda grid

`002_3_dem_table.py`
- make table of demographics

**Modeling/Analysis**
`003_condition_specific_networks.py`
- run RCFGL (`runRCFGL.py`) on each both conditions
- save precision matrices

`004_explore_networks.py`
- make adjacency matrices (uses functions from `analysis/utils/myDstream_functions.py`)
    - use tau threshold to filter out weak edges
- analyze networks non-visually
    - connectivity, etc.

`005_analyze_networks.py`
- find top largest components for each condition
- save metabolites to text files

`005_1_plot_networks.py`
- plot different graphs comparing between conditions
    - shared nodes/edges
    - unique nodes/edges

`005_2_component_pairs.py`
- look at shared graph components

`005_3_sphingo_sizes.py`
- look for patterns in sphingolipid sizes

`006_diff_enrich.R`
- enrichment analysis for 
    - shared/unique nodes
    - component nodes
    - shared/unique nodes for similar components

`007_hub_metab_correlation.R`
- make heatmap of how component metabolites correlate




