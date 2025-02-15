# squeue -l|grep cpu
# squeue -u wenhaoz1
# squeue|grep cpu|wc -l
# scontrol hold               # qhold
# scontrol release        # qrls
# scancel                          # qdel
# salloc    # Interactive job

# ll -h

cluster_head="wenhaoz1@mind.cs.cmu.edu:"
cluster_dir="/home/wenhaoz1/pjyu_temp_V1/SpkSims"
local_dir="/media/DATA1/Study/CompNeuro/Projects/Micro-clustering/Model_SigCorr/SpkSims"


ssh -Y wenhaoz1@mind.cs.cmu.edu
cd $cluster_dir
# mkdir Parameters Results
rm *.txt *.out

cd /media/DATA1/Study/CompNeuro/Projects/Functions_simul
scp histogram_mean_sem.m triu_new.m dtheta.m $cluster_head$cluster_dir
cd /media/DATA1/Study/CompNeuro/Projects/Micro-clustering/functions
scp FR_Analysis.m FR_Analysis_NoFitting.m $cluster_head$cluster_dir
cd $local_dir
scp EIF1DRFfastslowSyn_flexK.mexa64 spktime2count.mexa64 FFWD_spiking_generation.m $cluster_head$cluster_dir
scp Results_Processing.m Results_Processing.sh $cluster_head$cluster_dir
scp main.m main.sh $cluster_head$cluster_dir
cd $local_dir"/Parameters"
scp Parameters_FFWD.mat $cluster_head$cluster_dir"/Parameters"
scp Parameters_Recurrent_2.mat $cluster_head$cluster_dir"/Parameters"

# sbatch --array=1-20 main.sh
# sbatch --array=21-40 main.sh
# 
sbatch --array=2 Results_Processing.sh

scp $cluster_head$cluster_dir"/Results_Total_2.mat" $local_dir"/Results"





cd $local_dir"/Parameters"
scp Parameters_Recurrent_3.mat Parameters_Recurrent_4.mat Parameters_Recurrent_5.mat Parameters_Recurrent_6.mat Parameters_Recurrent_7.mat Parameters_Recurrent_8.mat Parameters_Recurrent_9.mat Parameters_Recurrent_10.mat Parameters_Recurrent_11.mat $cluster_head$cluster_dir"/Parameters"
#
sbatch --array=41-220%20 main.sh
# 
sbatch --array=3-11 Results_Processing.sh
#
zip Results_Total_1_11.zip Results_Total_*.mat
#
scp $cluster_head$cluster_dir"/Results_Total_1_11.zip" $local_dir"/Results"







cd $local_dir"/Parameters"
scp Parameters_Recurrent_12.mat Parameters_Recurrent_13.mat Parameters_Recurrent_14.mat Parameters_Recurrent_15.mat Parameters_Recurrent_16.mat Parameters_Recurrent_17.mat Parameters_Recurrent_18.mat Parameters_Recurrent_19.mat Parameters_Recurrent_20.mat $cluster_head$cluster_dir"/Parameters"
#
sbatch --array=221-400%21 main.sh
# 
sbatch --array=12-20 Results_Processing.sh
#
zip Results_Total_1_20.zip Results_Total_*.mat
#
scp $cluster_head$cluster_dir"/Results_Total_1_20.zip" $local_dir"/Results"




################################################################################3


cd $local_dir"/Parameters"
scp Parameters_Recurrent_6.mat Parameters_Recurrent_7.mat Parameters_Recurrent_8.mat Parameters_Recurrent_9.mat Parameters_Recurrent_10.mat Parameters_Recurrent_11.mat $cluster_head$cluster_dir"/Parameters"
#
sbatch --array=101-220%20 main.sh
# 
sbatch --array=6-11 Results_Processing.sh
#
zip Results_Total_6_11.zip Results_Total_6.mat Results_Total_7.mat Results_Total_8.mat Results_Total_9.mat Results_Total_10.mat Results_Total_11.mat
#
scp $cluster_head$cluster_dir"/Results_Total_6_11.zip" $local_dir"/Results"



cd $local_dir"/Parameters"
scp Parameters_Recurrent_15.mat Parameters_Recurrent_16.mat Parameters_Recurrent_17.mat Parameters_Recurrent_18.mat Parameters_Recurrent_19.mat Parameters_Recurrent_20.mat $cluster_head$cluster_dir"/Parameters"
#
sbatch --array=281-400%20 main.sh
# 
sbatch --array=15-20 Results_Processing.sh
#
zip Results_Total_15_20.zip Results_Total_15.mat Results_Total_16.mat Results_Total_17.mat Results_Total_18.mat Results_Total_19.mat Results_Total_20.mat
#
scp $cluster_head$cluster_dir"/Results_Total_15_20.zip" $local_dir"/Results"



