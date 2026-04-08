#!/bin/bash

# root=/Users/canderson/Documents/school/local-kechris-lab/kechris-lab/ # local root
ROOT=/projects/canderson2@xsede.org/kechris-lab/ # alpine root

cd "$ROOT/smoking-networks/" || exit 1

### Create Environment ###
conda env create -f environment.yaml

### Activate Environment ###
source "$(conda info --base)/etc/profile.d/conda.sh"  # workaround `conda init` for scripted initialization (i think)
conda activate smoknet-env

### Install other dependencies (mostly needed because of the oldish RCFGL build)###
# python -m pip install cppyy # installs to wrong location 
conda install -c conda-forge cppyy 

python -m pip install pycparser cffi # used for prox_tv build
# pip install  --no-build-isolation prox_tv
pip install prox_tv

### Install RFGL repo ###
dir="RCFGL"
[ -d "$dir" ]  && rm -rf "$dir" && echo "removed $dir"

git clone https://github.com/sealx017/RCFGL.git

echo -e "\n*** *** *** ***\\nRemember to change the ElasticNet normalize argument to manually normalizing with sklearn.preprocessing.Normalize\nlook here: \`analysis/utils/fixed_get_screening_2.py\`\n*** *** *** ***\n"

### OPTIONAL ###
### Install gnu parallel if needed and on mac ###
[ "$(uname -s)" == "Darwin" ] && brew install parallel


