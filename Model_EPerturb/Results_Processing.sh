#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=24GB
#SBATCH --time=0-00:09:00
#SBATCH --mail-user=yupeijia.qbio@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=Perturb_Sum

# Go to working directory
cd /home/wenhaoz1/pjyu_temp_V1/Perturb

module load matlab-9.11
matlab -nodisplay -nosplash -r Results_Processing 1>"Results_Processing_"$SLURM_ARRAY_TASK_ID"_out.txt" 2>"Results_Processing_"$SLURM_ARRAY_TASK_ID"_err.txt"
