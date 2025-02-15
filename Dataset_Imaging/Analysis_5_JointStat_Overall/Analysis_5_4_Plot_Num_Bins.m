if ispc, dirD = 'D:'; elseif isunix, dirD = '/media/DATA1'; end
addpath(genpath([dirD, '/Study/CompNeuro/Projects/Functions_simul/']));
dir0 = [dirD, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging'];
load([dir0, '/General_information.mat']);


d_BinStep_individual = 7.5;
for d_BinStep_overlap = [2.5 5]

dir_save = [dir0, '/Analysis_5_JointStat_Overall/JointStat/d_Bin', num2str(d_BinStep_overlap)];
%
load([dir0, '/Analysis_3_JointStat_EachDataset/Results_Datasets_dBin_',...
    num2str(d_BinStep_individual, '%.1f'), '_', num2str(d_BinStep_overlap, '%.1f'), '.mat']);
load([dir0, '/Analysis_5_JointStat_Overall/Results_Overall_dBin_', num2str(d_BinStep_individual,...
    '%.1f'), '_', num2str(d_BinStep_overlap, '%.1f'), '.mat']);
neuron_idx_dataset = [zeros(1, N_neuropilFactor); cumsum(N_neuron_OS_2_datasets, 1)];
%
if d_BinStep_overlap == 5
    idx_start = find(abs(d_BinCenter_overlap - 7.5) < 1e-3);
elseif d_BinStep_overlap == 2.5
    idx_start = find(abs(d_BinCenter_overlap - 7.5) < 1e-3);
end
%
Dataset_type_title = {'L2/3\_cyt.', 'L2/3\_nuc.', 'L4\_cyt.', 'L4\_som.'};
Dataset_type_title_full = {'L2/3\_cytosolic\_GCaMP', 'L2/3\_nuclear\_GCaMP', 'L4\_cytosolic\_GCaMP', 'L4\_soma\_GCaMP'};
Dataset_name_title = {'KF19', 'KF20', 'D12', 'D3', 'D4', 'Y18\_x0y0', 'Y18\_x-256y55',...
    'Y22\_x0y0', 'Y22\_x-401y418', 'Y24\_x0y0', 'Y24\_x426y133', 'Y25\_x0y0', 'Y25\_x10y-256', 'Y26', 'M199'};
%
hsv_clb = hsv * 0.875; mksize = 10; %scatter_edge_width = 0.25;
theta1 = 0; theta2 = 180; dtheta = 7.5; theta_edge = (theta1 - dtheta / 2): dtheta: (theta2 + dtheta / 2);
lwdth = 1.5; mksize = 12.5; cpsize = 7.5; txtsz = 15;
Clb_4 = [1 0.5 0; 1 0 0; 0 0 1; 0.5 0 1];


% %% Overlap results of each dataset
% %
% for i = 1: length(Dataset_type)
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.55, 1]); hold on;
% suptitle({['Number of neuron pairs per distance bin, ', Dataset_type_title{i}],...
        % ['(Bin Size = ', num2str(d_BinStep_overlap), ' \mum)']}, 8, 0.95);
% l = zeros(1, Dataset_N(i)); lgdtxt = cell(1, Dataset_N(i));
% clb = jet(Dataset_N(i)) * 0.875; clb = clb(end: -1: 1, :);
% for k = 1: Dataset_N(i)
    % dataset_k = Dataset_idx_edge(i) + k;
    % factor_k = neuropilFactor_best_idx(dataset_k);
    % %
    % if ROI_size_Tot(dataset_k, 1) < 300    % 256
        % [~, idx_end] = min(abs(d_BinCenter_overlap - 240));
    % else
        % [~, idx_end] = min(abs(d_BinCenter_overlap - 500));
    % end
    % %
    % x_data = d_BinCenter_overlap(idx_start: idx_end);
    % y_data = Num_Bin_datasets(idx_start: idx_end, dataset_k, factor_k);
    % y_data_se = SigCorr_se_datasets(idx_start: idx_end, dataset_k, factor_k);
    % idx_nonempty = find(y_data_se > eps);    % = 0, i.e. empty. This is possible after idx_start.
    % if length(idx_nonempty) ~= length(x_data)
        % x_data = x_data(idx_nonempty); y_data = y_data(idx_nonempty); y_data_se = y_data_se(idx_nonempty);
    % end
    % if Dataset_ToUse_SigCorr(dataset_k) == 1, lst = '-'; else, lst = '--'; end
    % l(k) = errorbar(x_data, y_data, (d_BinStep_overlap / 2) * ones(size(x_data)), 'horizontal', 'Marker', 'o',...
        % 'Color', clb(k, :), 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
    % lgdtxt{k} = Dataset_name_title{dataset_k};
% end
% set(gca, 'YScale', 'log'); axis square; grid on;
% axis([0 50 1 10000]); set(gca, 'XTick', 0: 5: 50);
% legend(l, lgdtxt, 'FontSize', txtsz);
% xlabel('Horizontal Cortical Distance (\mum)');
% set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
% %
% pause(2); print(gcf, '-dpng', [dir_save, '/Num_Bins_', Dataset_type{i}, '.png']);
% close;
% end


%% Overall Num_Bin of L2/3c, L2/3n, L4c
if d_BinStep_overlap == 5
    idx_start = find(abs(d_BinCenter_overlap - 7.5) < 1e-3);
elseif d_BinStep_overlap == 2.5
    idx_start = find(abs(d_BinCenter_overlap - 7.5) < 1e-3);
end
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.55, 1]); hold on;
suptitle({'Number of neuron pairs per distance bin',...
        ['(Bin Size = ', num2str(d_BinStep_overlap), ' \mum)']}, 8, 0.95);
l = zeros(1, length(Dataset_N)); lgdtxt = cell(1, length(Dataset_N));
for i = 1: length(Dataset_type)
    if (i == 1) | (i == 4)
        [~, idx_end] = min(abs(d_BinCenter_overlap - 500));
    else
        [~, idx_end] = min(abs(d_BinCenter_overlap - 240));
    end
    %
    x_data = d_BinCenter_overlap(idx_start: idx_end);
    y_data = Num_Bin_overall(idx_start: idx_end, i);
    y_data_se = SigCorr_se_overall(idx_start: idx_end, i);
    idx_nonempty = find(y_data_se > eps);    % = 0, i.e. empty. This is possible after idx_start.
    if length(idx_nonempty) ~= length(x_data)
        x_data = x_data(idx_nonempty); y_data = y_data(idx_nonempty);
    end
    l(i) = errorbar(x_data, y_data, (d_BinStep_overlap / 2) * ones(size(x_data)), 'horizontal', 'Marker', 'o',...
        'Color', Clb_4(i, :), 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
    lgdtxt{i} = Dataset_type_title_full{i};
end
set(gca, 'YScale', 'log'); axis square; grid on;
axis([0 50 1 10000]); set(gca, 'XTick', 0: 5: 50);
legend(l, lgdtxt, 'FontSize', txtsz, 'Location', 'southeast');
xlabel('Horizontal Cortical Distance (\mum)');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
pause(2); print(gcf, '-dpng', [dir_save, '/Num_Bins_Overall.png']);
close;


end





