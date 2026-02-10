#!/bin/zsh
##!/usr/bin/bash


### Go to Project Directory ###
# set working directory based on system
if [[ $OSTYPE == linux* ]]; then
  echo "OSTYPE Linux"
  os=linux
  dir=/projects/canderson2@xsede.org/kechris-lab/smoking-networks/analysis-versions/version001/
elif [[ $OSTYPE == darwin* ]]; then
  echo "OSTYPE Darwin"
  os=darwin
  dir=/Users/canderson/Documents/school/local-kechris-lab/kechris-lab/smoking-networks/analysis-versions/version001/
else 
  echo "OS Not Identified"
  exit 1
fi

cd $dir || exit 1
pwd 


# ### Activate Conda ###
# source "$HOME/miniconda3/etc/profile.d/conda.sh" # initialize for local
# # module purge # for alpine
# # module load anaconda  # for alpine
# conda activate smoknet-env

### Pipeline ###
# >>> Preprocess Data
F=analysis/001_preprocess_data.R
[[ -f $F ]] && echo -e "\n•••Running $F•••\n" || {echo -e "\nError: $F not found\n"; exit 1;}
conda run -n smoknet-env Rscript --vanilla $F
echo -e "\n•••$F Done•••\n"

# >>> Separate Data by smoking status
F=analysis/002_separate_conditions.py
[[ -f $F ]] && echo -e "\n•••Running $F•••\n" || {echo -e "\nError: $F not found\n"; exit 1;}
/Users/canderson/miniconda3/envs/smoknet-env/bin/python $F 
echo -e "\n•••$F Done•••\n"

# # >>> Run RCFGL lambda grid search
# F=analysis/002.1.1_slurm.sh
# echo Running $F
# echo -e "(Run \`squeeze\` to watch slurm job)"
# sbatch $F || exit 1
# echo -e "\n•••$F Done•••\n"

# # >>> Run lambda grid aic analysis
# F=analysis/002.2_analyze_aics.py
# echo Running $F
# bash $F || exit 1
# echo -e "\n•••$F Done•••\n"

# >>> Run RCFGL 
F=analysis/003_condition_specific_networks.py
[[ -f $F ]] && echo -e "\n•••Running $F•••\n" || {echo -e "\nError: $F not found\n"; exit 1;}
/Users/canderson/miniconda3/envs/smoknet-env/bin/python $F 
echo -e "\n•••$F Done•••\n"

# >>>> Analyze networks 1:
F=analysis.003_1_analyze_networks
echo -e "\n•••Running $F as module •••\n" 
/Users/canderson/miniconda3/envs/smoknet-env/bin/python -m "$F"
###### >>>> Process images
F=analysis/utils/process_images.sh
[[ -f $F ]] && echo -e "\n•••Running $F•••\n" || {echo -e "\nError: $F not found\n"; exit 1;}
zsh $F 
echo -e "\n•••$F Done•••\n"

# >>> >>> >>> >>> >>> >>> >>> >>> >>> >>> >>> >>> >>> >>> >>> >>> >>> 
echo -e "\n\\\\\\\\\\\\\n••• Pipeline Complete •••\n\\\\\\\\\\\\\n"
