#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8GB
#SBATCH --time=0-00:59:00
#SBATCH --mail-user=yupeijia.qbio@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=SigCorr_Sum

# Go to working directory
cd /home/wenhaoz1/pjyu_temp_V1/SpkSims

module load matlab-9.11
matlab -nodisplay -nosplash -r Results_Processing 1>"Results_Processing_"$SLURM_ARRAY_TASK_ID"_out.txt" 2>"Results_Processing_"$SLURM_ARRAY_TASK_ID"_err.txt"
