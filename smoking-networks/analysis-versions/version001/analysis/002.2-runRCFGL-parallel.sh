#!/bin/zsh

# number of processes to spawn
JOBS=20

# root=/Users/canderson/Documents/school/local-kechris-lab/kechris-lab/
ROOT=/projects/canderson2@xsede.org/kechris-lab/

cd "$ROOT/smoking-networks/analysis-versions/version001/" || exit 1

### Clean up result directory ###
[ -d results/002/aics ] && rm -rf 'results/002/aics/' && mkdir results/002/aics && echo 'removed results/002/aics'

## Activate conda env ###
# source "$(conda info --base)/etc/profile.d/conda.sh" 
source "$HOME/miniconda3/etc/profile.d/conda.sh" # initialize
conda activate smoknet-env # activate


### Make search grid ###
python << 'EOF'
import numpy as np
from itertools import product

vals = np.logspace(-3, 0, 20)
# vals = np.linspace(0,1, 20)

with open("results/002/lambda-grid.txt", "w") as f:
  for l1, l2 in product(vals, vals):
    f.write(f"{l1:.6g},{l2:.6g}\n")
EOF

### Run RCFGL in parallel and save to `aics` directory ###
# parallel -j 4 \
#   -j 4 \
#   python  analysis/runRCFGL.py --l1 {1} --l2 {2} \
#   ::: 0.001 0.01 0.1 1 \
#   ::: 0.001 0.01 0.1 1
  
parallel \
  -j $JOBS \
  --colsep ',' \
  /Users/canderson/miniconda3/envs/smoknet-env/bin/python  analysis/runRCFGL.py --l1 {1} --l2 {2} \
  :::: results/002/lambda-grid.txt
