#!/bin/bash

# root=/Users/canderson/Documents/school/local-kechris-lab/kechris-lab/
ROOT=/projects/canderson2@xsede.org/kechris-lab/

cd "$ROOT/smoking-networks/" || exit 1

### Create Environment ###
conda env create -f environment.yaml

### Activate Environment ###
source "$(conda info --base)/etc/profile.d/conda.sh"  # workaround `conda init`
conda activate smoknet-env

### Install other dependencies ###
# python -m pip install cppyy # installs to wrong location 
conda install -c conda-forge cppyy

python -m pip install pycparser cffi # used for prox_tv build
# pip install  --no-build-isolation prox_tv
pip install prox_tv

### Install RFGL repo ###
dir="RCFGL"
[ -d "$dir" ]  && rm -rf "$dir" && echo "removed $dir"

git clone https://github.com/sealx017/RCFGL.git

### Install gnu parallel if needed and on mac ###
[ "$(uname -s)" == "Darwin" ] && brew install parallel

### Install R dependencies ###
Rscript << 'EOF'
install.packages("BiocManager")
BiocManager::install('SummarizedExperiment')
install.packages("qs")
EOF