# squeue -l|grep cpu
# squeue -u wenhaoz1
# squeue|grep cpu|wc -l
# scontrol hold               # qhold
# scontrol release        # qrls
# scancel                          # qdel
# salloc    # Interactive job

# ll -h

ssh -Y wenhaoz1@mind.cs.cmu.edu
# cd $cluster_dir
# mkdir Parameters Results
rm *.txt *.out



cluster_head="wenhaoz1@mind.cs.cmu.edu:"
cluster_dir="/home/wenhaoz1/pjyu_temp_V1/NoiseCorr"
local_dir1="/media/DATA1/Study/CompNeuro/Projects/Micro-clustering/Model_NoiseCorr"
local_dir2="/media/DATA1/Study/CompNeuro/Projects/Micro-clustering/Model_SigCorr/SpkSims"

cd /media/DATA1/Study/CompNeuro/Projects/Functions_simul
scp histogram_mean_sem.m triu_new.m dtheta.m $cluster_head$cluster_dir
#
cd /media/DATA1/Study/CompNeuro/Projects/Micro-clustering/functions
scp FR_Analysis.m $cluster_head$cluster_dir
#
cd $local_dir1"/functions"
scp NoiseCorr_Analysis.m NoiseCorr_Analysis_dtheta.m $cluster_head$cluster_dir
#
cd $local_dir2
scp EIF1DRFfastslowSyn_flexK.mexa64 spktime2count.mexa64 FFWD_spiking_generation.m $cluster_head$cluster_dir
#
cd $local_dir1
scp FFWD_spiking_generation_LIF.m $cluster_head$cluster_dir
scp Results_Processing.m Results_Processing.sh $cluster_head$cluster_dir
scp main.m main.sh $cluster_head$cluster_dir
#
cd $local_dir1"/Parameters"
scp Parameters_FFWD.mat $cluster_head$cluster_dir"/Parameters"
scp Parameters_Recurrent_3.mat $cluster_head$cluster_dir"/Parameters"

sbatch --array=1-20,61-80%20 main.sh
sbatch --array=21-40,81-100%20 main.sh
sbatch --array=41-60,101-120%20 main.sh

sbatch --array=3,6 Results_Processing.sh

scp $cluster_head$cluster_dir"/Results_Total_3.mat" $local_dir1"/Results"


# scp $cluster_head$cluster_dir"/Results/Results_Par1_Simul1_tmp.mat" $local_dir1"/Results"


