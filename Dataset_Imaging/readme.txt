1. General information not included in raw data matrices, with General_information.m
2. Transfer raw data, in /Analysis_1_Individual_z
3. Find out all local pairs (d <= 15um) and manually pick them, in /Analysis_2_Local_Pairs
(done with automatic alpha)
4. Choose Neuropil alpha & use this individual dataset or not:
Correlation (joint statistics) of each dataset, and all alpha = [0, 0.25, 0.5, 0.75, 1, 1.25] and auto, in /Analysis_3_JointStat_EachDataset
Then Plot in /Analysis_4_NeuropilFactor_EachDataset
And determine which alpha to use for each dataset / use this dataset for total output or not, based on /Analysis_4_NeuropilFactor_EachDataset/d_Bin5/SigCorr, and in /Analysis_5_JointStat_Overall/Analysis_5_1_JointStat_Overall.m

