#!/bin/bash

#SBATCH --partition=amilan
#SBATCH --ntasks=10
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --qos=normal

#SBATCH --job-name=lambdaGS
#SBATCH --output=/projects/canderson2@xsede.org/kechris-lab/smoking-networks/analysis-versions/version001/logs/lambdaGS%j.out
#SBATCH --error=/projects/canderson2@xsede.org/kechris-lab/smoking-networks/analysis-versions/version001/logs/lambdaGS%j.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=christian.anderson@cuanschutz.edu


# --- Set working directory ---
dir=/projects/canderson2@xsede.org/kechris-lab/smoking-networks/analysis-versions/version001/
cd "$dir" || exit 1

# --- Make Logs directory ---
mkdir -p logs

# --- Load Anaconda and activate environment ---
module purge
module load anaconda
source "$(conda info --base)/etc/profile.d/conda.sh" 
conda activate smoknet-env

# --- Run script ---
SCRIPT="$dir/analysis/002.1-runRCFGL-parallel.sh"
bash $SCRIPT --jobs $SLURM_NTASKS
