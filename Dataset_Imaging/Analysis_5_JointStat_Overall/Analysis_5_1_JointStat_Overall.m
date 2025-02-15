if ispc, dirD = 'D:'; elseif isunix, dirD = '/media/DATA1'; end
addpath(genpath([dirD, '/Study/CompNeuro/Projects/Functions_simul/']));
dir0 = [dirD, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging'];
addpath([dir0, '/Analysis_5_JointStat_Overall']);
load([dir0, '/General_information.mat']);
load([dir0, '/Analysis_2_Local_Pairs/idx_remove.mat'], 'idx_remove');

%%%%%%%%%%%%%%%%
% How to integrate datasets into type overall -- Based on figure outputs in /Analysis_4_NeuropilFactor_EachDataset
neuropilFactor_best_idx = [[5, 5], [7, 7, 7], [5, 5, 5, 5, 5, 5, 5, 5, 5], [6]];    % 0, 0.25, 0.5, 0.75, 1, 1.25, auto
% Dataset_ToUse_SigCorr = [[0, 1], [1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 1], [1]];
% Dataset_ToUse_dtheta = [[1, 1], [1, 0, 1], [1, 1, 1, 1, 1, 1, 1, 1, 1], [1]];
Dataset_ToUse_SigCorr = [[1, 1], [1, 1, 1], [1, 1, 1, 1, 1, 1, 1, 1, 1], [1]];
Dataset_ToUse_dtheta = [[1, 1], [1, 0, 1], [1, 1, 1, 1, 1, 1, 1, 1, 1], [1]];
%%%%%%%%%%%%%%%%

d_BinStep_individual = 7.5;
for d_BinStep_overlap = [2.5 5]
% Building up overlap joint stat.
load([dir0, '/Analysis_3_JointStat_EachDataset/Results_Datasets_dBin_', num2str(d_BinStep_individual,...
    '%.1f'), '_', num2str(d_BinStep_overlap, '%.1f'), '.mat'], 'N_Bin', 'd_BinCenter_overlap', 'd_BinEdge_overlap');
%
SigCorr_mean_layer = NaN(N_Bin, 2); SigCorr_se_layer = NaN(N_Bin, 2);    % Total L2/3 or L4.
ROI_centeroid_OS_1_layer = []; dFF0_mean_OS_layer = []; Block_N_OS_1_layer = [];
%
SigCorr_mean_overall = NaN(N_Bin, length(Dataset_type));
SigCorr_se_overall = NaN(N_Bin, length(Dataset_type));
dtheta_pref_mean_overall = NaN(N_Bin, length(Dataset_type));
dtheta_pref_se_overall = NaN(N_Bin, length(Dataset_type));
Num_Bin_overall = NaN(N_Bin, length(Dataset_type));
%
for i = 1: length(Dataset_type)
ROI_centeroid_OS_1_o = []; dFF0_mean_OS_o = []; Block_N_OS_1_o = [];
ROI_centeroid_OS_2_o = []; Pref_orientation_OS_o = []; Block_N_OS_2_o = [];
%
for dataset_k = (Dataset_idx_edge(i) + 1): Dataset_idx_edge(i + 1)
z_value_list = Dataset_z{dataset_k};
factor_k = neuropilFactor_best_idx(dataset_k);
for z_k = 1: length(z_value_list)
    savename = [Dataset_type{i}, '_', Dataset_name{dataset_k}, '_z', num2str(z_value_list(z_k))];
    load([dir0, '/Analysis_1_Individual_z/', savename, '_s2p_NpSize30_Ana1.mat'],...
        'N_neuron', 'ROI_centeroid', 'dFF0_mean', 'Pref_orientation', 'isOS');
    % remove "invalid" neurons
    idx_remove_k = idx_remove{dataset_k}{z_k};    % Comes from neuropil factor auto. No difference...
    idx_valid = 1: N_neuron; idx_valid(idx_remove_k) = NaN; idx_valid = idx_valid(~isnan(idx_valid));
    ROI_centeroid = ROI_centeroid(idx_valid, :);
    dFF0_mean = dFF0_mean(idx_valid, :, factor_k);
    Pref_orientation = Pref_orientation(idx_valid, factor_k);
    isOS = isOS(idx_valid, :, factor_k);
    % pick up OS
    idxOS_1 = find(isOS(:, 1) == 1); idxOS_2 = find(isOS(:, 2) == 1);
    ROI_centeroid_OS_1 = ROI_centeroid(idxOS_1, :); dFF0_mean_OS = dFF0_mean(idxOS_1, :);
    ROI_centeroid_OS_2 = ROI_centeroid(idxOS_2, :); Pref_orientation_OS = Pref_orientation(idxOS_2);
    % recording
    if Dataset_ToUse_SigCorr(dataset_k) == 1
        ROI_centeroid_OS_1_o = [ROI_centeroid_OS_1_o; ROI_centeroid_OS_1];
        dFF0_mean_OS_o = [dFF0_mean_OS_o; dFF0_mean_OS];
        Block_N_OS_1_o = [Block_N_OS_1_o, length(idxOS_1)];
    end
    if Dataset_ToUse_dtheta(dataset_k) == 1
        ROI_centeroid_OS_2_o = [ROI_centeroid_OS_2_o; ROI_centeroid_OS_2];
        Pref_orientation_OS_o = [Pref_orientation_OS_o; Pref_orientation_OS];
        Block_N_OS_2_o = [Block_N_OS_2_o, length(idxOS_2)];
    end
end
end
%
[SigCorr_mean_overall(:, i), SigCorr_se_overall(:, i), Num_Bin_overall(:, i)] =...
    SignalCorr_Dist(dFF0_mean_OS_o, ROI_centeroid_OS_1_o, Block_N_OS_1_o, d_BinEdge_overlap);
[dtheta_pref_mean_overall(:, i), dtheta_pref_se_overall(:, i), ~] =...
    dtheta_pref_Dist(Pref_orientation_OS_o, ROI_centeroid_OS_2_o, Block_N_OS_2_o, d_BinEdge_overlap);
%
ROI_centeroid_OS_1_layer = [ROI_centeroid_OS_1_layer; ROI_centeroid_OS_1_o];
dFF0_mean_OS_layer = [dFF0_mean_OS_layer; dFF0_mean_OS_o];
Block_N_OS_1_layer = [Block_N_OS_1_layer, Block_N_OS_1_o];
if mod(i, 2) == 0    % i == 2 or 4
    [SigCorr_1d, dist_1d] = SigCorr_Integration(dFF0_mean_OS_layer, ROI_centeroid_OS_1_layer, Block_N_OS_1_layer);
    [SigCorr_mean_layer(:, i/2), SigCorr_se_layer(:, i/2), ~] = histogram_mean_sem(SigCorr_1d, dist_1d, d_BinEdge_overlap);
    ROI_centeroid_OS_1_layer = []; dFF0_mean_OS_layer = []; Block_N_OS_1_layer = [];
    clear SigCorr_1d dist_1d
end
end
%
save([dir0, '/Analysis_5_JointStat_Overall/Results_Overall_dBin_', num2str(d_BinStep_individual,...
    '%.1f'), '_', num2str(d_BinStep_overlap, '%.1f'), '.mat'], 'neuropilFactor_best_idx',...
    'Dataset_ToUse_SigCorr', 'Dataset_ToUse_dtheta', 'N_Bin',...
    'd_BinCenter_overlap', 'd_BinEdge_overlap', 'SigCorr_mean_overall', 'SigCorr_se_overall',...
    'dtheta_pref_mean_overall', 'dtheta_pref_se_overall', 'Num_Bin_overall',...
    'SigCorr_mean_layer', 'SigCorr_se_layer');
end

