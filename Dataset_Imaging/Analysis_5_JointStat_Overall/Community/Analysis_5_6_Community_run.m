if ispc, dirD = 'D:'; elseif isunix, dirD = '/media/DATA1'; end
addpath(genpath([dirD, '/Study/CompNeuro/Projects/Functions_simul/']));
dir0 = [dirD, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging'];

addpath([dir0, '/Analysis_5_JointStat_Overall/Community']);
load([dir0, '/General_information.mat']);
load([dir0, '/Analysis_2_Local_Pairs/idx_remove.mat'], 'idx_remove');
load([dir0, '/Analysis_5_JointStat_Overall/Results_Overall_dBin_7.5_5.0.mat'],...
    'neuropilFactor_best_idx', 'Dataset_ToUse_SigCorr');

dist_thr = 15;    % um
dtheta_thr = 30;    % deg
N_trial_shuffle = 50;

ROI_centeroid_OS_Tot = cell(1, length(Dataset_type));
Pref_orientation_OS_Tot = cell(1, length(Dataset_type));
Block_N_OS_Tot = cell(1, length(Dataset_type));
N_OS_Tot = zeros(1, length(Dataset_type));
%
for i = 1: length(Dataset_type)
ROI_centeroid_OS_2_o = []; Pref_orientation_OS_o = []; Block_N_OS_2_o = [];
%
for dataset_k = (Dataset_idx_edge(i) + 1): Dataset_idx_edge(i + 1)
if Dataset_ToUse_SigCorr(dataset_k) == 1
    z_value_list = Dataset_z{dataset_k};
    factor_k = neuropilFactor_best_idx(dataset_k);
    for z_k = 1: length(z_value_list)
        savename = [Dataset_type{i}, '_', Dataset_name{dataset_k}, '_z', num2str(z_value_list(z_k))];
        load([dir0, '/Analysis_1_Individual_z/', savename, '_s2p_NpSize30_Ana1.mat'],...
            'N_neuron', 'ROI_centeroid', 'Pref_orientation', 'isOS');
        % remove "invalid" neurons
        idx_remove_k = idx_remove{dataset_k}{z_k};    % Comes from neuropil factor auto. No difference...
        idx_valid = 1: N_neuron; idx_valid(idx_remove_k) = NaN; idx_valid = idx_valid(~isnan(idx_valid));
        ROI_centeroid = ROI_centeroid(idx_valid, :);
        Pref_orientation = Pref_orientation(idx_valid, factor_k);
        isOS = isOS(idx_valid, :, factor_k);
        % pick up OS
        idxOS_2 = find(isOS(:, 2) == 1);
        ROI_centeroid_OS_2 = ROI_centeroid(idxOS_2, :);
        Pref_orientation_OS = Pref_orientation(idxOS_2);
        % recording
        ROI_centeroid_OS_2_o = [ROI_centeroid_OS_2_o; ROI_centeroid_OS_2];
        Pref_orientation_OS_o = [Pref_orientation_OS_o; Pref_orientation_OS];
        Block_N_OS_2_o = [Block_N_OS_2_o, length(idxOS_2)];
    end
end
end
ROI_centeroid_OS_Tot{i} = ROI_centeroid_OS_2_o;
Pref_orientation_OS_Tot{i} = Pref_orientation_OS_o;
Block_N_OS_Tot{i} = Block_N_OS_2_o;
N_OS_Tot(i) = length(Pref_orientation_OS_o);
end
clear ROI_centeroid_OS_2 Pref_orientation_OS idxOS_2
clear ROI_centeroid_OS_2_o Pref_orientation_OS_o Block_N_OS_2_o
%N_OS_Tot    % [2014, 911, 2441, 1124]
% Experimental conditions are the same for L2/3 nuc. and L4 cyt., but not for L/3 cyt. or L4 soma.
% You can choose only L2/3 nuc. and L4 cyt., so they will have similar ratio_within_comm.


for i = 1: length(Dataset_type)
    ROI_centeroid_OS_o = ROI_centeroid_OS_Tot{i};
    Pref_orientation_OS_o = Pref_orientation_OS_Tot{i};
    Block_N_OS_o = Block_N_OS_Tot{i};
    %
    Comm_result(i).ROI_centeroid_OS = ROI_centeroid_OS_o;
    Comm_result(i).Pref_orientation_OS = Pref_orientation_OS_o;
    [~, Comm_result(i).Comm_neuron_idx, Comm_result(i).Size_comm,...
        Comm_result(i).idx_within_comm, Comm_result(i).idx_beyond_comm] = Find_clusters_func(...
        ROI_centeroid_OS_o, Pref_orientation_OS_o, Block_N_OS_o, dist_thr, dtheta_thr, 0);
    %
    Size_comm_shuffled = cell(1, N_trial_shuffle);
    N_within_comm_shuffled = zeros(1, N_trial_shuffle);
    N_beyond_comm_shuffled = zeros(1, N_trial_shuffle);
    trial_k = 0;
    while trial_k < N_trial_shuffle
    try
        trial_k = trial_k + 1;
        [Pref_orientation_OS_o_shuffled_eg, Comm_neuron_idx_shuffled_eg, Size_comm_shuffled{trial_k},...
            idx_within_comm_shuffled_eg, idx_beyond_comm_shuffled_eg] = Find_clusters_func(...
            ROI_centeroid_OS_o, Pref_orientation_OS_o, Block_N_OS_o, dist_thr, dtheta_thr, 1);
        N_within_comm_shuffled(trial_k) = length(idx_within_comm_shuffled_eg);
        N_beyond_comm_shuffled(trial_k) = length(idx_beyond_comm_shuffled_eg);
        fprintf([num2str(trial_k), ' / ', num2str(N_trial_shuffle), ' fin.\n']);
    catch
        trial_k = trial_k - 1; fprintf('oops...\n'); continue;
    end
    end
    ratio_within_comm_shuffled = N_within_comm_shuffled ./...
        (N_within_comm_shuffled + N_beyond_comm_shuffled);
    %
    % save
    Comm_result(i).Pref_orientation_OS_o_shuffled_eg = Pref_orientation_OS_o_shuffled_eg;
    Comm_result(i).Comm_neuron_idx_shuffled_eg = Comm_neuron_idx_shuffled_eg;
    Comm_result(i).Size_comm_shuffled = Size_comm_shuffled;
    Comm_result(i).idx_within_comm_shuffled_eg = idx_within_comm_shuffled_eg;
    Comm_result(i).idx_beyond_comm_shuffled_eg = idx_beyond_comm_shuffled_eg;
    Comm_result(i).ratio_within_comm_shuffled = ratio_within_comm_shuffled;
    %
    fprintf([Dataset_type{i}, ' fin.\n']);
end


savename = [dir0, '/Analysis_5_JointStat_Overall/Community/Community_distThr_',...
        num2str(dist_thr), '_dthetaThr_', num2str(dtheta_thr), '.mat'];
save(savename, 'N_OS_Tot', 'dist_thr', 'dtheta_thr', 'N_trial_shuffle', 'Comm_result');




for i = 1: length(Dataset_type)
    a = length(Comm_result(i).idx_within_comm);
    b = length(Comm_result(i).idx_beyond_comm);
    r1 = a / (a + b);
    r2 = mean(Comm_result(i).ratio_within_comm_shuffled);
    d2 = std(Comm_result(i).ratio_within_comm_shuffled);
    fprintf([Dataset_type{i}, ', raw: ', num2str(r1, '%.2f'), ' (', num2str(a),...
        ' / ', num2str(a + b), '), shuffled: ', num2str(r2, '%.2f'), ' +/- ', num2str(d2, '%.2f'), '.\n']);
end
