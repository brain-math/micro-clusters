if ispc, dirD = 'D:'; elseif isunix, dirD = '/media/DATA1'; end
dir0 = [dirD, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging'];
dir_name = [dir0, '/Analysis_1_Individual_z']; if ~exist(dir_name, 'dir'), mkdir(dir_name); end
dir_name = [dir0, '/Analysis_2_Local_Pairs']; if ~exist(dir_name, 'dir'), mkdir(dir_name); end
dir_name = [dir0, '/Analysis_3_JointStat_EachDataset']; if ~exist(dir_name, 'dir'), mkdir(dir_name); end
dir_name = [dir0, '/Analysis_4_NeuropilFactor_EachDataset']; if ~exist(dir_name, 'dir'), mkdir(dir_name); end
dir_name = [dir0, '/Analysis_5_JointStat_Overall']; if ~exist(dir_name, 'dir'), mkdir(dir_name); end

Dataset_type = {'L23_cytosolic_GCaMP', 'L23_nuclear_GCaMP', 'L4_cytosolic_GCaMP', 'L4_soma_GCaMP'};
Dataset_N = [2, 3, 9, 1]; Dataset_N_tot = sum(Dataset_N);
Dataset_idx_edge = [0 cumsum(Dataset_N)];
Dataset_name = {'KF19', 'KF20',...
    'D12', 'D3', 'D4',...
    'Y18_x0y0', 'Y18_x-256y55', 'Y22_x0y0', 'Y22_x-401y418', 'Y24_x0y0', 'Y24_x426y133', 'Y25_x0y0', 'Y25_x10y-256', 'Y26',...
    'M199'};
%
Dataset_z = cell(1, Dataset_N_tot);
N_frame = cell(1, Dataset_N_tot);
% KF19
Dataset_z{1} = [175, 200, 225];
N_frame{1} = repmat([35, 5, 30], length(Dataset_z{1}), 1);    % [Total frame #, pre. interstm. frame #, evoked frame #]
% KF20
Dataset_z{2} = [114, 134, 154, 174, 200, 225, 250];
N_frame{2} = repmat([35, 5, 30], length(Dataset_z{2}), 1);
% D12   
Dataset_z{3} = [180, 200, 220, 234, 250, 270, 290];
N_frame{3} = [repmat([16, 3, 10], 4, 1); repmat([13, 3, 10], 3, 1)];
% D3
Dataset_z{4} = [120, 140, 160, 180, 200, 220];
N_frame{4} = repmat([24, 6, 18], length(Dataset_z{4}), 1);
% D4
Dataset_z{5} = [120, 140, 160, 180, 200];
N_frame{5} = repmat([24, 6, 18], length(Dataset_z{5}), 1);
% Y18_x0y0
Dataset_z{6} = [320, 340, 360, 380, 400, 420, 440];
N_frame{6} = repmat([24, 6, 18], length(Dataset_z{6}), 1);
% Y18_x-256y55
Dataset_z{7} = [340, 360, 380, 400, 420, 440, 460];
N_frame{7} = repmat([24, 6, 18], length(Dataset_z{7}), 1);
% Y22_x0y0
Dataset_z{8} = [320, 340, 360, 380, 400];
N_frame{8} = repmat([24, 6, 18], length(Dataset_z{8}), 1);
% Y22_x-401y418
Dataset_z{9} = [260, 280, 300, 320, 340, 360, 380];
N_frame{9} = [repmat([22, 5, 17], 3, 1); [24, 6, 18]; repmat([22, 5, 17], 3, 1)];
% Y24_x0y0
Dataset_z{10} = [260, 280, 300, 320, 340];
N_frame{10} = repmat([22, 5, 17], length(Dataset_z{10}), 1);
% Y24_x426y133
Dataset_z{11} = [300, 320, 340, 360];
N_frame{11} = repmat([22, 5, 17], length(Dataset_z{11}), 1);
% Y25_x0y0
Dataset_z{12} = [280, 300, 320];
N_frame{12} = repmat([22, 5, 17], length(Dataset_z{12}), 1);
% Y25_x10y-256
Dataset_z{13} = [220, 240, 260];
N_frame{13} = repmat([22, 5, 17], length(Dataset_z{13}), 1);
% Y26
Dataset_z{14} = [320, 340, 360];
N_frame{14} = repmat([24, 6, 18], length(Dataset_z{14}), 1);
% M199
Dataset_z{15} = [210, 230, 250, 270, 290, 310, 330, 350];
N_frame{15} = repmat([30, 5, 25], length(Dataset_z{15}), 1);
% [Total frame #, pre. interstm. frame #, evoked frame #]


xls_title = {'ROI index', 'fittedAngle', 'pvalue', 'DSI_fitted', 'DSI_raw',...
    'OSI_fitted', 'OSI_raw', 'gDSI', 'gOSI', 'theta(pref_direction)',...
    'FWHM', 'maxY', 'minY_avgForFit', 'minOptions(0:N/A;1:shift;2:negSet0)', 'isOS?',...
    'OS_maxDf_f>10%', 'OS_pvalue<0.05, used', 'OS_goodFitting (SSE<0.40, RSQ>0.60)', 'fitted OSI>0.00', 'isDS?',...
    'DS_maxDf_f>10%', 'DS_pvalue<0.05, used', 'DS_goodFittingSSE<0.40RSQ>0.60', 'fittedDSI>0.50', 'DS_ isOS? used'};

% The frame rate, L2/3 cyt., L2/3 nuc., L4 cyt., L4 som., unit: frame / s
frame_rate = [7.5, 2.18, 2.18, 7.5];

% 12 angles, each 10 trials
N_trial = 10; theta_stim = linspace(0, 330, 12); N_theta = length(theta_stim);

% micron / pixel
pixel2um = [0.821, 0.821,...
    1, 1, 1,...
    1, 1, 1, 1, 1, 1, 1, 1, 1,...
    1];

% Standard of "OS" neuron, we have two:
    % (i) Non-flatten (Pvalue_anova < 0.05) & Responsive (max DFF0 > 0.1)
    % (ii) Non-flatten & Responsive & OS_goodFitting (SSE < 0.40, RSQ > 0.60) (i.e. must be well-shaped two-peak curve)
    % (i) used for signal correlation; (ii) used for anything related to pref. orientation.
Thr_Pvalue_anova = 0.05; Thr_Responsive = 0.1;

neuropilFactor_manual = [0, 0.25, 0.5, 0.75, 1, 1.25];
N_neuropilFactor = length(neuropilFactor_manual) + 1;

save('General_information.mat', 'Dataset_type', 'Dataset_N', 'Dataset_N_tot', 'Dataset_idx_edge',...
    'Dataset_name', 'Dataset_z', 'N_frame', 'xls_title', 'frame_rate', 'N_trial', 'theta_stim', 'N_theta',...
    'pixel2um', 'Thr_Pvalue_anova', 'Thr_Responsive', 'neuropilFactor_manual', 'N_neuropilFactor');

