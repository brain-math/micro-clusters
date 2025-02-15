if ispc, dirD = 'D:'; elseif isunix, dirD = '/media/DATA1'; end
dir0 = [dirD, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging'];

% KF19: 25
% KF20: only the 2nd part, 25
% D12: cut into two parts (delete 234), 20
% D3: 20
% Do not use D4.
% All L4: 20
Dataset_type = {'L23_cytosolic_GCaMP', 'L23_nuclear_GCaMP', 'L4_cytosolic_GCaMP', 'L4_soma_GCaMP'};
Dataset_N = [2, 3, 9, 1]; Dataset_N_tot = sum(Dataset_N);
Dataset_idx_edge = [0 cumsum(Dataset_N)];
Dataset_name = {'KF19', 'KF20', 'D12', 'D12', 'D3', 'Y18_x0y0', 'Y18_x-256y55',...
    'Y22_x0y0', 'Y22_x-401y418', 'Y24_x0y0', 'Y24_x426y133', 'Y25_x0y0', 'Y25_x10y-256', 'Y26', 'M199'};
%
Dataset_z = cell(1, Dataset_N_tot);
N_frame = cell(1, Dataset_N_tot);
%
Dataset_z{1} = [175, 200, 225];    % KF19
Dataset_z{2} = [200, 225, 250];    % KF20, 2
Dataset_z{3} = [180, 200, 220];    % D12, 1
Dataset_z{4} = [250, 270, 290];    % D12, 2
Dataset_z{5} = [120, 140, 160, 180, 200, 220];    % D3
Dataset_z{6} = [320, 340, 360, 380, 400, 420, 440];    % Y18_x0y0
Dataset_z{7} = [340, 360, 380, 400, 420, 440, 460];    % Y18_x-256y55
Dataset_z{8} = [320, 340, 360, 380, 400];    % Y22_x0y0
Dataset_z{9} = [260, 280, 300, 320, 340, 360, 380];    % Y22_x-401y418
Dataset_z{10} = [260, 280, 300, 320, 340];    % Y24_x0y0
Dataset_z{11} = [300, 320, 340, 360];    % Y24_x426y133
Dataset_z{12} = [280, 300, 320];    % Y25_x0y0
Dataset_z{13} = [220, 240, 260];    % Y25_x10y-256
Dataset_z{14} = [320, 340, 360];    % Y26
Dataset_z{15} = [210, 230, 250, 270, 290, 310, 330, 350];    % M199

theta_stim = linspace(0, 330, 12); N_theta = length(theta_stim);

load([dir0, '/Analysis_5_JointStat_Overall/Results_Overall_dBin_7.5_5.0.mat'], 'neuropilFactor_best_idx');


load([dir0, '/Analysis_2_Local_Pairs/idx_remove.mat'], 'idx_remove');
idx_remove_old = idx_remove; clear idx_remove
idx_remove = cell(1, Dataset_N_tot);
idx_remove{1} = idx_remove_old{1};    % KF19
idx_remove{2} = idx_remove_old{2}(5: 7);    % KF20, 2
idx_remove{3} = idx_remove_old{3}(1: 3);    % D12, 1
idx_remove{4} = idx_remove_old{3}(5: 7);    % D12, 2
idx_remove{5} = idx_remove_old{4};    % D3
for i = 6: Dataset_N_tot, idx_remove{i} = idx_remove_old{i}; end; clear i idx_remove_old


save([dir0, '/Analysis_5_JointStat_Overall/Column/General_information_Z.mat'],...
    'Dataset_idx_edge', 'Dataset_N', 'Dataset_N_tot', 'Dataset_name', 'Dataset_type',...
    'Dataset_z', 'idx_remove', 'N_frame', 'N_theta', 'neuropilFactor_best_idx', 'theta_stim');


