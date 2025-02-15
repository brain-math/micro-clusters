#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8GB
#SBATCH --time=0-03:29:59
#SBATCH --mail-user=yupeijia.qbio@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name=Perturb

# Go to working directory
cd /home/wenhaoz1/pjyu_temp_V1/Perturb

# module avail
module load matlab-9.11
matlab -nodisplay -nosplash -r main 1>"main_"$SLURM_ARRAY_TASK_ID"_out.txt" 2>"main_"$SLURM_ARRAY_TASK_ID"_err.txt"
