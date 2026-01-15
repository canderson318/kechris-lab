#!/bin/bash

### Create Environment ###
conda env create -f environment.yaml

### Activate Environment ###
# workaround `conda init`
source "$(conda info --base)/etc/profile.d/conda.sh" 
conda activate smoknet-env

### Install other dependencies ###
python -m pip install cppyy
python -m pip install pycparser cffi # used for prox_tv build
# pip install  --no-build-isolation prox_tv
pip install prox_tv


