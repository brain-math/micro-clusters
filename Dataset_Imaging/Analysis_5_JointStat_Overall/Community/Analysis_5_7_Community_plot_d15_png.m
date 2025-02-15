if ispc, dirD = 'D:'; elseif isunix, dirD = '/media/DATA1'; end
addpath(genpath([dirD, '/Study/CompNeuro/Projects/Functions_simul/']));
dir0 = [dirD, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging'];
load([dir0, '/General_information.mat']);
load([dir0, '/Analysis_3_JointStat_EachDataset/Results_Datasets_dBin_7.5_5.0.mat'], 'ROI_size_Tot');
ROI_size_Tot = ROI_size_Tot(Dataset_idx_edge(1: end - 1) + 1, :);

dist_thr = 15;    % um
dtheta_thr = 30;    % deg
%
load([dir0, '/Analysis_5_JointStat_Overall/Community/Community_distThr_',...
    num2str(dist_thr), '_dthetaThr_', num2str(dtheta_thr), '.mat'],...
    'N_OS_Tot', 'N_trial_shuffle', 'Comm_result');

Dataset_type_title = {'L2/3\_cytosolic\_GCaMP', 'L2/3\_nuclear\_GCaMP',...
    'L4\_cytosolic\_GCaMP', 'L4\_soma\_GCaMP'};
%
hsv_clb = hsv * 0.875; mksize = 7.5; mksize_gray = 2.5;
%mksize = 10; mksize_gray = 3;
scatter_edge_clr = 0.6 * ones(1, 3); line_width = [0.25, 0.75, 0.75];
%
x_ctr = 2: 3; x_edge = [1.5: 1: 3.5];


figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.7, 1]);
suptitle(['Distance threshold = ', num2str(dist_thr), ' \mum, \Delta\theta\_pref. threshold = ',...
    num2str(dtheta_thr), ' deg.'], 4, 0.98);
%
for i = 1: length(Dataset_type)
    ROI_centeroid_OS = Comm_result(i).ROI_centeroid_OS;
    Pref_orientation_OS = Comm_result(i).Pref_orientation_OS;
    Comm_neuron_idx = Comm_result(i).Comm_neuron_idx;
    Size_comm = Comm_result(i).Size_comm;
    idx_within_comm = Comm_result(i).idx_within_comm;
    idx_beyond_comm = Comm_result(i).idx_beyond_comm;
    %
    Pref_orientation_OS_shuffled_eg = Comm_result(i).Pref_orientation_OS_o_shuffled_eg;
    Comm_neuron_idx_shuffled_eg = Comm_result(i).Comm_neuron_idx_shuffled_eg;
    Size_comm_shuffled = Comm_result(i).Size_comm_shuffled;
    idx_within_comm_shuffled_eg = Comm_result(i).idx_within_comm_shuffled_eg;
    idx_beyond_comm_shuffled_eg = Comm_result(i).idx_beyond_comm_shuffled_eg;
    ratio_within_comm_shuffled = Comm_result(i).ratio_within_comm_shuffled;
    %
    ratio_within_comm = length(idx_within_comm) / (length(idx_within_comm) + length(idx_beyond_comm));
    ratio_within_comm_shuffled_mean = mean(ratio_within_comm_shuffled);
    ratio_within_comm_shuffled_std = std(ratio_within_comm_shuffled);
    %
    histc_raw = histcounts(Size_comm, x_edge, 'normalization', 'count');
    histc_shuffle = zeros(N_trial_shuffle, length(x_ctr));
    histc_ratio_raw2shuffle = zeros(N_trial_shuffle, length(x_ctr));
    for k = 1: N_trial_shuffle
        histc_shuffle(k, :) = histcounts(Size_comm_shuffled{k}, x_edge, 'normalization', 'count');
        histc_ratio_raw2shuffle(k, :) = histc_raw ./ histc_shuffle(k, :);
    end
    histc_ratio_raw2shuffle(abs(histc_ratio_raw2shuffle) > 1145141919810) = NaN;
    histc_ratio_raw2shuffle_mean = mean(histc_ratio_raw2shuffle, 1, 'omitnan');
    histc_ratio_raw2shuffle_std = std(histc_ratio_raw2shuffle, [], 1, 'omitnan');


    N_neuron_k = length(Pref_orientation_OS);
    color_idx = ceil((Pref_orientation_OS / 180) * size(hsv_clb, 1)); color_idx(color_idx == 0) = 1;
    color_array = zeros(N_neuron_k, 3);
    for k = 1: N_neuron_k, color_array(k, :) = hsv_clb(color_idx(k), :); end
    ROI_centeroid_OS_x = ROI_centeroid_OS(:, 2);
    ROI_centeroid_OS_y = ROI_size_Tot(i, 2) - ROI_centeroid_OS(:, 1);
    %
    subplot(3, length(Dataset_type), i); hold on;
    scatter(ROI_centeroid_OS_x(idx_beyond_comm), ROI_centeroid_OS_y(idx_beyond_comm),...
        mksize_gray, 'MarkerEdgeColor', scatter_edge_clr);
    scatter(ROI_centeroid_OS_x(idx_within_comm), ROI_centeroid_OS_y(idx_within_comm), mksize,...
        color_array(idx_within_comm, :), 'fill', 'MarkerEdgeColor', 'none');
    for k = 1: length(Comm_neuron_idx)
        idx_tmp = Comm_neuron_idx{k};
        for ii = 1: (length(idx_tmp) - 1)
        for jj = (ii + 1): length(idx_tmp)
            if length(idx_tmp) == 2, width_k = 1; elseif length(idx_tmp) == 3, width_k = 2; else, width_k = 3; end
            plot([ROI_centeroid_OS_x(idx_tmp(ii)), ROI_centeroid_OS_x(idx_tmp(jj))],...
                [ROI_centeroid_OS_y(idx_tmp(ii)), ROI_centeroid_OS_y(idx_tmp(jj))],...
                'Color', 'k', 'LineWidth', line_width(width_k));
        end
        end
    end
    axis square; axis([0, ROI_size_Tot(i, 1), 0, ROI_size_Tot(i, 2)]);
    set(gca, 'XTick', linspace(0, ROI_size_Tot(i, 1), 9), 'YTick', linspace(0, ROI_size_Tot(i, 2), 9));
    title({Dataset_type_title{i}, [num2str(ratio_within_comm * 100, '%.1f'), '% (',...
        num2str(length(idx_within_comm)), ' / ', num2str(N_neuron_k), ')']}, 'FontWeight', 'normal');
    set(gca, 'box', 'on', 'TickDir', 'out');


    color_idx = ceil((Pref_orientation_OS_shuffled_eg / 180) * size(hsv_clb, 1)); color_idx(color_idx == 0) = 1;
    color_array = zeros(N_neuron_k, 3);
    for k = 1: N_neuron_k, color_array(k, :) = hsv_clb(color_idx(k), :); end
    %
    subplot(3, length(Dataset_type), i + length(Dataset_type)); hold on;
    scatter(ROI_centeroid_OS_x(idx_beyond_comm_shuffled_eg),...
        ROI_centeroid_OS_y(idx_beyond_comm_shuffled_eg),...
        mksize_gray, 'MarkerEdgeColor', scatter_edge_clr);
    scatter(ROI_centeroid_OS_x(idx_within_comm_shuffled_eg),...
        ROI_centeroid_OS_y(idx_within_comm_shuffled_eg), mksize,...
        color_array(idx_within_comm_shuffled_eg, :), 'fill', 'MarkerEdgeColor', 'none');
    for k = 1: length(Comm_neuron_idx_shuffled_eg)
        idx_tmp = Comm_neuron_idx_shuffled_eg{k};
        for ii = 1: (length(idx_tmp) - 1)
        for jj = (ii + 1): length(idx_tmp)
            if length(idx_tmp) == 2, width_k = 1; elseif length(idx_tmp) == 3, width_k = 2; else, width_k = 3; end
            plot([ROI_centeroid_OS_x(idx_tmp(ii)), ROI_centeroid_OS_x(idx_tmp(jj))],...
                [ROI_centeroid_OS_y(idx_tmp(ii)), ROI_centeroid_OS_y(idx_tmp(jj))],...
                'Color', 'k', 'LineWidth', line_width(width_k));
        end
        end
    end
    axis square; axis([0, ROI_size_Tot(i, 1), 0, ROI_size_Tot(i, 2)]);
    set(gca, 'XTick', linspace(0, ROI_size_Tot(i, 1), 9), 'YTick', linspace(0, ROI_size_Tot(i, 2), 9));
    title({'Shuffled', [num2str(ratio_within_comm_shuffled_mean * 100, '%.1f'), '% (std: ',...
        num2str(ratio_within_comm_shuffled_std * 100, '%.1f'), '%)']}, 'FontWeight', 'normal');
    set(gca, 'box', 'on', 'TickDir', 'out');
    

    subplot(3, length(Dataset_type), i + length(Dataset_type) * 2); hold on;
    bar(x_ctr, histc_ratio_raw2shuffle_mean, 0.75, 'EdgeColor', 'k', 'FaceColor', scatter_edge_clr);
    errorbar(x_ctr, histc_ratio_raw2shuffle_mean, histc_ratio_raw2shuffle_std, 'Color', 'k', 'LineStyle', 'none');
    for k = 1: length(x_ctr)
        text(x_ctr(k), histc_ratio_raw2shuffle_mean(k) + histc_ratio_raw2shuffle_std(k),...
            num2str(histc_raw(k)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
    end
    plot([0 6], [1 1], 'k--');
    axis square; axis([min(x_edge) - 0.3, max(x_edge) + 0.3, 0, 7.1]);
    set(gca, 'XTick', x_ctr, 'YTick', 0: 7);
    if i == 1, xlabel('Number of neurons in clusters');
    title('Counts, real / Counts, shuffled', 'FontWeight', 'normal'); end
end


savename = ['Community_distThr_', num2str(dist_thr), '_dthetaThr_', num2str(dtheta_thr)];
pause(2); print(gcf, '-dpng', [dir0, '/Analysis_5_JointStat_Overall/Community/', savename, '.png']);
close;

