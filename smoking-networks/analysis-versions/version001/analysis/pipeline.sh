#!/usr/bin/bash

### Go to Project Directory ###
dir='/projects/canderson2@xsede.org/kechris-lab/smoking-networks/analysis-versions/version001/'
cd $dir || exit 1

### Activate Conda ###
# source "$HOME/miniconda3/etc/profile.d/conda.sh" # initialize\
module purge
module load anaconda
conda activate smoknet-env

### Pipeline ###
# >>> Preprocess Data
F=analysis/001-preprocess-data.R
echo -e "\n•••Running $F•••\n"
Rscript --vanilla $F # > /dev/null
echo -e "\n•••$F Done•••\n"

# >>> Separate Data by smoking status
F=analysis/002-separate-conditions.py
echo -e "\n•••Running $F•••\n"
python $F 
echo -e "\n•••$F Done•••\n"

# >>> Run RCFGL lambda grid search
# F=analysis/002.1.1-slurm.sh
# echo Running $F
# sbatch $F
# echo -e "\n•••$F Done•••\n"

# # >>> Run lambda grid aic analysis
# F=analysis/002.2-analyze-aics.py
# echo Running $F
# bash $F
# echo -e "\n•••$F Done•••\n"


echo -e "\n\\\\\\\\\\\\\n••• Pipeline Complete •••\n\\\\\\\\\\\\\n"

echo -e "Run \`squeeze\` to watch slurm job"