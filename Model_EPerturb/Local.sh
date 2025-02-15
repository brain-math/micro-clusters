# squeue -l|grep cpu
# squeue -u wenhaoz1
# squeue|grep cpu|wc -l
# scontrol hold               # qhold
# scontrol release        # qrls
# scancel                          # qdel
# salloc    # Interactive job

# ll -h

cluster_head="wenhaoz1@mind.cs.cmu.edu:"
cluster_dir="/home/wenhaoz1/pjyu_temp_V1/Perturb"
local_dir="/media/DATA1/Study/CompNeuro/Projects/Micro-clustering/Model_EPerturb"


ssh -Y wenhaoz1@mind.cs.cmu.edu
# mkdir Parameters Results

cd $local_dir
scp EIF1DRFfastslowSyn_flexK.mexa64 EIF1DRFfastslowSyn_flexK_PertE.mexa64 spktime2count.mexa64 FFWD_spiking_generation.m $cluster_head$cluster_dir
scp Results_Processing.m Results_Processing.sh $cluster_head$cluster_dir
scp main.m main.sh $cluster_head$cluster_dir
cd $local_dir"/Parameters"
scp Parameters_Recurrent_4.mat $cluster_head$cluster_dir"/Parameters"

# sbatch --array=1-301%38 main.sh
# sbatch --array=302-602%22 main.sh
# sbatch --array=603-903%22 main.sh
sbatch --array=904-1504%31 main.sh
# 
sbatch --array=4 Results_Processing.sh

scp $cluster_head$cluster_dir"/Results_Total_4_A.mat" $local_dir"/Results"
scp $cluster_head$cluster_dir"/Results_Total_4_B.mat" $local_dir"/Results"

# ls -1v|cat -n > results.log

