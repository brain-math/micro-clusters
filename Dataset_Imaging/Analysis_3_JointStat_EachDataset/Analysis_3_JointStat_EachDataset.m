if ispc, dirD = 'D:'; elseif isunix, dirD = '/media/DATA1'; end
addpath(genpath([dirD, '/Study/CompNeuro/Projects/Functions_simul/']));
dir0 = [dirD, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging'];

load([dir0, '/General_information.mat']);
load([dir0, '/Analysis_2_Local_Pairs/idx_remove.mat'], 'idx_remove');

%% INTEGRATE ALL .MAT FILES.
d_BinStep_individual = 7.5;
for d_BinStep_overlap = [2.5 5]

N_z = 0; for k = 1: Dataset_N_tot, N_z = N_z + length(Dataset_z{k}); end    % 72
%
% For preferred orientation maps:
ROI_size_Tot = NaN(Dataset_N_tot, 2);
N_neuron_OS_2_individual = NaN(N_z, N_neuropilFactor);    % Each z.
N_neuron_overall_2_individual = NaN(N_z, N_neuropilFactor);
N_neuron_OS_2_datasets = NaN(Dataset_N_tot, N_neuropilFactor);    % Each dataset.
N_neuron_overall_2_datasets = NaN(Dataset_N_tot, N_neuropilFactor);
% idx = [0 cumsum(N)]; (idx(k) + 1: idx(k + 1))
ROI_centeroid_OS_2_Tot = cell(1, N_neuropilFactor);    % Each (N_neuron_Tot_OS_2, 2) in order.
Pref_orientation_OS_Tot = cell(1, N_neuropilFactor);    % Each (N_neuron_Tot_OS_2, 1) in order.
% 
% For joint stat:
N_Bin = 201;
d_BinCenter_individual = [[0: 1: (N_Bin - 2)]' * d_BinStep_individual;...
    ((N_Bin - 0.5) * d_BinStep_individual + 1500) / 2];
d_BinEdge_individual = [([0: 1: (N_Bin - 1)] - 0.5) * d_BinStep_individual, 1500];
%
if d_BinStep_overlap == 2.5
    d_BinCenter_overlap = [[0: 1: (N_Bin - 2)]' * d_BinStep_overlap;...
        ((N_Bin - 0.5) * d_BinStep_overlap + 1500) / 2];
    d_BinEdge_overlap = [([0: 1: (N_Bin - 1)] - 0.5) * d_BinStep_overlap, 1500];
elseif d_BinStep_overlap == 5
    d_BinCenter_overlap = [([0: 1: (N_Bin - 2)]' + 0.5) * d_BinStep_overlap; 1250];
    d_BinEdge_overlap = [[0: 1: (N_Bin - 1)] * d_BinStep_overlap, 1500];
end
%
name = {'SigCorr_mean', 'SigCorr_se', 'dtheta_pref_mean', 'dtheta_pref_se', 'Num_Bin'};
% Num_Bin follows isOS(:, 1), since served for joint stat.
for k = 1: length(name)
    eval([name{k}, '_individual = NaN(N_Bin, N_z, N_neuropilFactor);']);
    eval([name{k}, '_datasets = NaN(N_Bin, Dataset_N_tot, N_neuropilFactor);']);
end; clear name k


for factor_k = 1: N_neuropilFactor
% Accumulating all types for each neuropil factor.
ROI_centeroid_OS_2_f = []; Pref_orientation_OS_f = [];
%
dataset_k = 0; z_j = 0;    % z_j here is total j among all z_depth .mat files.
for i = 1: length(Dataset_type)
    % Accumulating all dataset for each overall type, L2/3c, L2/3n, L4c.
    ROI_centeroid_OS_1_o = []; dFF0_mean_OS_o = [];
    ROI_centeroid_OS_2_o = []; Pref_orientation_OS_o = [];
    Block_N_OS_1_o = []; Block_N_OS_2_o = [];
    for j = 1: Dataset_N(i)
        dataset_k = dataset_k + 1;
        z_value_list = Dataset_z{dataset_k};
        % Accumulating all z for each dataset.
        ROI_centeroid_OS_1_d = []; dFF0_mean_OS_d = [];
        ROI_centeroid_OS_2_d = []; Pref_orientation_OS_d = [];
        Block_N_OS_1_d = zeros(1, length(z_value_list)); Block_N_OS_2_d = zeros(1, length(z_value_list));
        N_neuron_overall_2_datasets_tmp = 0;
        for z_k = 1: length(z_value_list)
            z_j = z_j + 1;
            savename = [Dataset_type{i}, '_', Dataset_name{dataset_k}, '_z', num2str(z_value_list(z_k))];
            load([dir0, '/Analysis_1_Individual_z/', savename, '_s2p_NpSize30_Ana1.mat'],...
                'N_neuron', 'ROI_size', 'ROI_centeroid', 'dFF0_mean', 'Pref_orientation', 'isOS');
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
            N_neuron_OS_2_individual(z_j, factor_k) = length(idxOS_2);
            N_neuron_overall_2_individual(z_j, factor_k) = length(idx_valid);
            N_neuron_overall_2_datasets_tmp = N_neuron_overall_2_datasets_tmp + length(idx_valid);
            ROI_centeroid_OS_1_d = [ROI_centeroid_OS_1_d; ROI_centeroid_OS_1];
            dFF0_mean_OS_d = [dFF0_mean_OS_d; dFF0_mean_OS];
            ROI_centeroid_OS_2_d = [ROI_centeroid_OS_2_d; ROI_centeroid_OS_2];
            Pref_orientation_OS_d = [Pref_orientation_OS_d; Pref_orientation_OS];
            Block_N_OS_1_d(z_k) = length(idxOS_1); Block_N_OS_2_d(z_k) = length(idxOS_2);
            % Joint stat for individual z
            [SigCorr_mean_individual(:, z_j, factor_k), SigCorr_se_individual(:, z_j, factor_k),...
                Num_Bin_individual(:, z_j, factor_k)] = SignalCorr_Dist(dFF0_mean_OS,...
                ROI_centeroid_OS_1, Block_N_OS_1_d(z_k), d_BinEdge_individual);
            [dtheta_pref_mean_individual(:, z_j, factor_k), dtheta_pref_se_individual(:, z_j, factor_k), ~] =...
                dtheta_pref_Dist(Pref_orientation_OS, ROI_centeroid_OS_2, Block_N_OS_2_d(z_k), d_BinEdge_individual);
            %
            fprintf(['z', num2str(z_value_list(z_k)), ' fin.    ']);
        end; fprintf('\n');
        % Joint stat for each dataset.
        [SigCorr_mean_datasets(:, dataset_k, factor_k), SigCorr_se_datasets(:, dataset_k, factor_k),...
            Num_Bin_datasets(:, dataset_k, factor_k)] = SignalCorr_Dist(dFF0_mean_OS_d,...
            ROI_centeroid_OS_1_d, Block_N_OS_1_d, d_BinEdge_overlap);
        [dtheta_pref_mean_datasets(:, dataset_k, factor_k), dtheta_pref_se_datasets(:, dataset_k, factor_k), ~] =...
            dtheta_pref_Dist(Pref_orientation_OS_d, ROI_centeroid_OS_2_d, Block_N_OS_2_d, d_BinEdge_overlap);
        % recording
        N_neuron_OS_2_datasets(dataset_k, factor_k) = size(ROI_centeroid_OS_2_d, 1);
        N_neuron_overall_2_datasets(dataset_k, factor_k) = N_neuron_overall_2_datasets_tmp;
        ROI_centeroid_OS_1_o = [ROI_centeroid_OS_1_o; ROI_centeroid_OS_1_d];
        dFF0_mean_OS_o = [dFF0_mean_OS_o; dFF0_mean_OS_d];
        ROI_centeroid_OS_2_o = [ROI_centeroid_OS_2_o; ROI_centeroid_OS_2_d];
        Pref_orientation_OS_o = [Pref_orientation_OS_o; Pref_orientation_OS_d];
        Block_N_OS_1_o = [Block_N_OS_1_o, Block_N_OS_1_d];
        Block_N_OS_2_o = [Block_N_OS_2_o, Block_N_OS_2_d];
        ROI_size_Tot(dataset_k, :) = ROI_size;
        %
        fprintf([Dataset_name{dataset_k}, ' fin.\n']);
    end
    %
    ROI_centeroid_OS_2_f = [ROI_centeroid_OS_2_f; ROI_centeroid_OS_2_o];
    Pref_orientation_OS_f = [Pref_orientation_OS_f; Pref_orientation_OS_o];
    fprintf([Dataset_type{i}, ' fin.\n']);
end
%
ROI_centeroid_OS_2_Tot{factor_k} = ROI_centeroid_OS_2_f;
Pref_orientation_OS_Tot{factor_k} = Pref_orientation_OS_f;
%
if factor_k ~= N_neuropilFactor
    fprintf(['neuropilFactor ', num2str(neuropilFactor_manual(factor_k)), ' fin.\n']);
else
    fprintf('neuropilFactor auto fin.\n');
end
end


neuropilFactor_auto_avg_individual = NaN(1, N_z);
neuropilFactor_auto_avg_datasets = zeros(1, Dataset_N_tot);
%
dataset_k = 0; z_j = 0;
for i = 1: length(Dataset_type)
for j = 1: Dataset_N(i)
dataset_k = dataset_k + 1;
z_value_list = Dataset_z{dataset_k};
for z_k = 1: length(z_value_list)
    z_j = z_j + 1;
    savename = [Dataset_type{i}, '_', Dataset_name{dataset_k}, '_z', num2str(z_value_list(z_k))];
    load([dir0, '/RawData/', savename, '_s2p_NpSize30.mat'], 'neuropilFactor_auto');
    neuropilFactor_auto_avg_individual(z_j) = mean(neuropilFactor_auto);
    % Actually identical for each z. So don't worry about idx_remove.
    neuropilFactor_auto_avg_datasets(dataset_k) = neuropilFactor_auto_avg_datasets(dataset_k) +...
        neuropilFactor_auto_avg_individual(z_j) / length(z_value_list);
end
end
end


save([dir0, '/Analysis_3_JointStat_EachDataset/Results_Datasets_dBin_',...
    num2str(d_BinStep_individual, '%.1f'), '_', num2str(d_BinStep_overlap, '%.1f'), '.mat'],...
    'ROI_size_Tot', 'N_neuron_OS_2_individual', 'N_neuron_overall_2_individual',...
    'N_neuron_OS_2_datasets', 'N_neuron_overall_2_datasets',...
    'ROI_centeroid_OS_2_Tot', 'Pref_orientation_OS_Tot', 'N_Bin',...
    'd_BinCenter_individual', 'd_BinCenter_overlap', 'd_BinEdge_individual', 'd_BinEdge_overlap',...
    'SigCorr_mean_individual', 'SigCorr_se_individual', 'dtheta_pref_mean_individual',...
    'dtheta_pref_se_individual', 'Num_Bin_individual',...
    'SigCorr_mean_datasets', 'SigCorr_se_datasets', 'dtheta_pref_mean_datasets',...
    'dtheta_pref_se_datasets', 'Num_Bin_datasets',...
    'neuropilFactor_auto_avg_individual', 'neuropilFactor_auto_avg_datasets');

end

