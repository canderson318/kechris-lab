#!/usr/bin/bash

# number of processes to spawn
JOBS=1   # default

while [[ $# -gt 0 ]]; do
    case "$1" in
        --jobs|-j)
            JOBS="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

echo "Using JOBS=$JOBS"

# root=/Users/canderson/Documents/school/local-kechris-lab/kechris-lab/
dir=/projects/canderson2@xsede.org/kechris-lab/

cd "$dir/smoking-networks/analysis-versions/version001/" || exit 1

### Clean up result directory ###
[ -d results/002/aics ] && rm -rf 'results/002/aics/' && mkdir results/002/aics && echo 'removed results/002/aics'

## Activate conda env ###
# source "$(conda info --base)/etc/profile.d/conda.sh" 
# conda activate smoknet-env # activate


### Make search grid ###
NUM_VALS=20
NUM_VALS=$NUM_VALS python << 'EOF'
import numpy as np
import os
from itertools import product

# vals = np.logspace(-3, 0, int(os.environ["NUM_VALS"]))
vals = np.linspace(0,1, int(os.environ["NUM_VALS"]))

with open("results/002/lambda-grid.txt", "w") as f:
  for l1, l2 in product(vals, vals):
    f.write(f"{l1:.6g},{l2:.6g}\n")
EOF

### Run RCFGL in parallel and save to `aics` directory ###
# load parallel
module load gnu_parallel

# run
echo Running...

parallel \
  -j $JOBS \
  --colsep ',' \
  python analysis/runRCFGL.py --l1 {1} --l2 {2} \
  :::: results/002/lambda-grid.txt

echo Done