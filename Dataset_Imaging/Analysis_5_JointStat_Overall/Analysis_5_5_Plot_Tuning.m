if ispc, dirD = 'D:'; elseif isunix, dirD = '/media/DATA1'; end
addpath(genpath([dirD, '/Study/CompNeuro/Projects/Functions_simul/']));
dir0 = [dirD, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging'];
load([dir0, '/General_information.mat']);
load([dir0, '/Analysis_2_Local_Pairs/idx_remove.mat'], 'idx_remove');
dir_save = [dir0, '/Analysis_5_JointStat_Overall/Tuning']; if ~exist(dir_save, 'dir'), mkdir(dir_save); end
%
load([dir0, '/Analysis_3_JointStat_EachDataset/Results_Datasets_dBin_7.5_5.0.mat']);
load([dir0, '/Analysis_5_JointStat_Overall/Results_Overall_dBin_7.5_5.0.mat']);
neuron_idx_dataset = [zeros(1, N_neuropilFactor); cumsum(N_neuron_OS_2_datasets, 1)];
%
Dataset_type_title = {'L2/3\_cyt.', 'L2/3\_nuc.', 'L4\_cyt.', 'L4\_som.'};
Dataset_type_title_full = {'L2/3\_cytosolic\_GCaMP', 'L2/3\_nuclear\_GCaMP', 'L4\_cytosolic\_GCaMP', 'L4\_soma\_GCaMP'};
Dataset_name_title = {'KF19', 'KF20', 'D12', 'D3', 'D4', 'Y18\_x0y0', 'Y18\_x-256y55',...
    'Y22\_x0y0', 'Y22\_x-401y418', 'Y24\_x0y0', 'Y24\_x426y133', 'Y25\_x0y0', 'Y25\_x10y-256', 'Y26', 'M199'};
%
lwdth = 1.5; mksize = 12.5; cpsize = 7.5; txtsz = 15;
hsv_clb = hsv * 0.875; mksize = 10; %scatter_edge_width = 0.25;
%
theta1 = 0; theta2 = 180; dtheta = 7.5; theta_edge = (theta1 - dtheta / 2): dtheta: (theta2 + dtheta / 2);
%
x1 = 0; x2 = 1; dx = 0.02; x_edge = (x1 - dx / 2): dx: (x2 + dx / 2);
ytickg = cell(1, 8);
ytickg{1} = [0: 100: 600]; ytickg{2} = [0: 25: 150]; ytickg{3} = [0: 40: 240]; ytickg{4} = [0: 200: 1200];
ytickg{5} = [0: 20: 120]; ytickg{6} = [0: 10: 60]; ytickg{7} = [0: 25: 150]; ytickg{8} = [0: 10: 60];


If_shrink_L4soma = 1;

%% Pref. orientation maps
ROI_centeroid_OS_f = cell(1, length(Dataset_type));
Pref_orientation_OS_f = cell(1, length(Dataset_type));
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
k = 0;
for i = 1: length(Dataset_type)
ROI_centeroid_OS_2_o = []; Pref_orientation_OS_o = [];
%
for dataset_k = (Dataset_idx_edge(i) + 1): Dataset_idx_edge(i + 1)
if Dataset_ToUse_SigCorr(dataset_k) == 1
    k = k + 1;
    factor_k = neuropilFactor_best_idx(dataset_k);
    idx_k = [neuron_idx_dataset(dataset_k, factor_k) + 1: neuron_idx_dataset(dataset_k + 1, factor_k)];
    %
    ROI_centeroid_k = ROI_centeroid_OS_2_Tot{factor_k}(idx_k, :);
    Pref_orientation_k = Pref_orientation_OS_Tot{factor_k}(idx_k);
    N_neuron_k = N_neuron_OS_2_datasets(dataset_k, factor_k);
    %
    ROI_centeroid_OS_2_o = [ROI_centeroid_OS_2_o; ROI_centeroid_k];
    Pref_orientation_OS_o = [Pref_orientation_OS_o; Pref_orientation_k];
    %
    color_idx = ceil((Pref_orientation_k / 180) * size(hsv_clb, 1)); color_idx(color_idx == 0) = 1;
    color_array = zeros(N_neuron_k, 3);
    for neuron_k = 1: N_neuron_k, color_array(neuron_k, :) = hsv_clb(color_idx(neuron_k), :); end
    %
    subplot(3, 5, k);
    scatter(ROI_centeroid_k(:, 2), ROI_size_Tot(dataset_k, 2) - ROI_centeroid_k(:, 1),...    % xy
        mksize, color_array, 'fill', 'MarkerEdgeColor', 'none');%, 'LineWidth', scatter_edge_width); 
    axis square;
    if ROI_size_Tot(dataset_k, 1) < 300    % 256
        axis([0, ROI_size_Tot(dataset_k, 1), 0, ROI_size_Tot(dataset_k, 2)]);
        set(gca, 'XTick', 0: 32: 256, 'YTick', 0: 32: 256);
    elseif ROI_size_Tot(dataset_k, 1) < 1000    % 1024 * 0.821, 840
        axis([0, ROI_size_Tot(dataset_k, 1), 0, ROI_size_Tot(dataset_k, 2)]);
        set(gca, 'XTick', 0: 100: 800, 'YTick', 0: 100: 800);
    else    % 1024
        if If_shrink_L4soma == 0
            axis([0, ROI_size_Tot(dataset_k, 1), 0, ROI_size_Tot(dataset_k, 2)]);
            set(gca, 'XTick', 0: 200: 1000, 'YTick', 0: 200: 1000);
        elseif If_shrink_L4soma == 1
            axis([0, 512, 256, 768]);
            set(gca, 'XTick', 0: 64: 512, 'YTick', 256: 64: 768);
        end
    end
    if k == 1, xlabel('x (\mum)'); ylabel('y (\mum)'); end
    if k == Dataset_N_tot
        colormap(hsv_clb); clb = colorbar; clb.Position = [0.75, 0.11, 0.01, 0.217];
        set(clb, 'YTick', 0: 1/4: 1, 'YTickLabel', {'0^o', '45^o', '90^o', '135^o', '180^o'});
    end
    t1 = title(''); text(t1.Position(1), t1.Position(2) * 1.05,...
        [Dataset_type_title{i}, ', ', Dataset_name_title{dataset_k}, ' (N = ', num2str(N_neuron_k), ')'],...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontName', 'Helvetica', 'FontSize', 11);
    set(gca, 'box', 'on', 'TickDir', 'out');
end
end
%
ROI_centeroid_OS_f{i} = ROI_centeroid_OS_2_o;
Pref_orientation_OS_f{i} = Pref_orientation_OS_o;
end
%
pause(2);
if If_shrink_L4soma == 0, print(gcf, '-dpng', [dir_save, '/Pref_orientation_Maps_a.png']);
elseif If_shrink_L4soma == 1, print(gcf, '-dpng', [dir_save, '/Pref_orientation_Maps_b.png']); end
close;
%
%
ROI_size_Tot_i = ROI_size_Tot(Dataset_idx_edge(1: end - 1) + 1, :);
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 0.5]);
for i = 1: length(Dataset_type)
    ROI_centeroid_k = ROI_centeroid_OS_f{i};
    N_neuron_k = size(ROI_centeroid_k, 1);
    %
    color_idx = ceil((Pref_orientation_OS_f{i} / 180) * size(hsv_clb, 1)); color_idx(color_idx == 0) = 1;
    color_array = zeros(N_neuron_k, 3);
    for neuron_k = 1: N_neuron_k, color_array(neuron_k, :) = hsv_clb(color_idx(neuron_k), :); end
    %
    subplot(1, length(Dataset_type), i);
    scatter(ROI_centeroid_k(:, 2), ROI_size_Tot_i(i, 2) - ROI_centeroid_k(:, 1),...    % xy
        mksize, color_array, 'fill', 'MarkerEdgeColor', 'none');%, 'LineWidth', scatter_edge_width); 
    axis square;
    if ROI_size_Tot_i(i, 1) < 300    % 256
        axis([0, ROI_size_Tot_i(i, 1), 0, ROI_size_Tot_i(i, 2)]);
        set(gca, 'XTick', 0: 32: 256, 'YTick', 0: 32: 256);
    elseif ROI_size_Tot_i(i, 1) < 1000    % 1024 * 0.821, 840
        axis([0, ROI_size_Tot_i(i, 1), 0, ROI_size_Tot_i(i, 2)]);
        set(gca, 'XTick', 0: 100: 800, 'YTick', 0: 100: 800);
    else    % 1024
        if If_shrink_L4soma == 0
            axis([0, ROI_size_Tot(dataset_k, 1), 0, ROI_size_Tot(dataset_k, 2)]);
            set(gca, 'XTick', 0: 200: 1000, 'YTick', 0: 200: 1000);
        elseif If_shrink_L4soma == 1
            axis([0, 512, 256, 768]);
            set(gca, 'XTick', 0: 64: 512, 'YTick', 256: 64: 768);
        end
    end
    if i == 1, xlabel('x (\mum)'); ylabel('y (\mum)'); end
    if i == length(Dataset_type)
        colormap(hsv_clb); clb = colorbar; clb.Position = [0.92, 0.18, 0.01, 0.675];
        set(clb, 'YTick', 0: 1/4: 1, 'YTickLabel', {'0^o', '45^o', '90^o', '135^o', '180^o'});
    end
    t1 = title(''); text(t1.Position(1), t1.Position(2) * 1.05, [Dataset_type_title{i}, ' (N = ', num2str(N_neuron_k), ')'],...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontName', 'Helvetica', 'FontSize', 11);
    set(gca, 'FontSize', txtsz, 'box', 'on', 'TickDir', 'out');
end
%
pause(2);
if If_shrink_L4soma == 0, print(gcf, '-dpng', [dir_save, '/Pref_orientation_Maps_Overall_a.png']);
elseif If_shrink_L4soma == 1, print(gcf, '-dpng', [dir_save, '/Pref_orientation_Maps_Overall_b.png']); end
close;



%% Distribution of Preferred Orientations
pref_theta_histcounts = zeros(length(Dataset_type), length(theta_edge) - 1);
title_txt = cell(1, length(Dataset_type));
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
k = 0;
for i = 1: length(Dataset_type)
Pref_orientation_k_tot = []; N_OS = 0; N_overall = 0;
for dataset_k = (Dataset_idx_edge(i) + 1): Dataset_idx_edge(i + 1)
if Dataset_ToUse_SigCorr(dataset_k) == 1
    k = k + 1;
    factor_k = neuropilFactor_best_idx(dataset_k);
    idx_k = [neuron_idx_dataset(dataset_k, factor_k) + 1: neuron_idx_dataset(dataset_k + 1, factor_k)];
    Pref_orientation_k = Pref_orientation_OS_Tot{factor_k}(idx_k);
    Pref_orientation_k_tot = [Pref_orientation_k_tot; Pref_orientation_k];
    N_neuron_k = N_neuron_OS_2_datasets(dataset_k, factor_k);
    N_OS = N_OS + N_neuron_k;
    N_overall = N_overall + N_neuron_overall_2_datasets(dataset_k, factor_k);
    %
    subplot(3, 5, k);
    histogram(Pref_orientation_k, theta_edge, 'Normalization', 'count', 'FaceColor', [0 0.5 1]);
    [pval, ~] = circ_otest(Pref_orientation_k * 2 * (pi / 180), dtheta * (pi / 180));
    if pval >= 0.01, title2 = ['p = ', num2str(pval, '%.3f')]; else, title2 = ['p < 10^{', num2str(ceil(log10(pval))), '}']; end
    axis square; xlim([theta1 theta2]); set(gca, 'XTick', 0: 15: 180, 'XTickLabel',...
        {'0', '', '30', '', '60', '', '90', '', '120', '', '150', '', '180'});
    title({[Dataset_type_title{i}, ', ', Dataset_name_title{dataset_k}],...
        [title2, ' (N = ', num2str(N_neuron_k), ')']}, 'FontWeight', 'normal');
    if k == 1
        xlabel('Pref. orientation'); ylabel('Number count');
    end
    if k == 2
        ylim([0 200]); set(gca, 'YTick', 0: 40: 200);
    elseif k == 14
        ylim([0 100]); set(gca, 'YTick', 0: 20: 100);
    else
        ylim([0 60]); set(gca, 'YTick', 0: 10: 60);
    end
    set(gca, 'box', 'off', 'TickDir', 'out');
end
end
%
y = histcounts(Pref_orientation_k_tot, theta_edge, 'Normalization', 'count');
y([1 end]) = ((y(1) + y(end)) / 2) / 2;
pref_theta_histcounts(i, :) = y;
%
[pval, ~] = circ_otest(Pref_orientation_k_tot * 2 * (pi / 180), dtheta * (pi / 180));
if pval >= 0.01
    title_txt{i} = ['N = ', num2str(N_OS), ' / ', num2str(N_overall), ', p = ', num2str(pval, '%.3f')];
else
    title_txt{i} = ['N = ', num2str(N_OS), ' / ', num2str(N_overall), ', p < 10^{', num2str(ceil(log10(pval))), '}'];
end
end
%
pause(2); print(gcf, '-dpng', [dir_save, '/Pref_orientation_Distribution.png']);
close;
%
ymax = 200; dy = 50;
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 0.6]);
for i = 1: length(Dataset_type)
    subplot(1, length(Dataset_type), i);
    bar(theta1: dtheta: theta2, pref_theta_histcounts(i, :), 1, 'FaceColor', [0 0.5 1]);
    axis square; axis([theta1 theta2 0 ymax]); set(gca, 'XTick', 0: 15: 180, 'XTickLabel',...
        {'0', '', '30', '', '60', '', '90', '', '120', '', '150', '', '180'}, 'YTick', 0: dy: ymax);
    title([Dataset_type_title{i}, ', ', title_txt{i}],  'FontWeight', 'normal');
    if i == 1
        xlabel('Pref. orientation'); ylabel('Number count');
    end
    set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
end
pause(2); print(gcf, '-dpng', [dir_save, '/Pref_orientation_Distribution_Overall.png']);
close;


%% gOSI
gOSI_Tot = cell(1, length(Dataset_type));
for i = 1: length(Dataset_type)
gOSI_Tot_k = []; gOSI_Tot_OS_k = [];
for dataset_k = (Dataset_idx_edge(i) + 1): Dataset_idx_edge(i + 1)
if Dataset_ToUse_SigCorr(dataset_k) == 1
    z_value_list = Dataset_z{dataset_k};
    factor_k = neuropilFactor_best_idx(dataset_k);
    for z_k = 1: length(z_value_list)
        savename = [Dataset_type{i}, '_', Dataset_name{dataset_k}, '_z', num2str(z_value_list(z_k))];
        load([dir0, '/RawData/', savename, '_s2p_NpSize30.mat'], 'gOSI');
        load([dir0, '/Analysis_1_Individual_z/', savename, '_s2p_NpSize30_Ana1.mat'], 'N_neuron', 'isOS');
        idx_remove_k = idx_remove{dataset_k}{z_k};    % Comes from neuropil factor auto. No difference...
        idx_valid = 1: N_neuron; idx_valid(idx_remove_k) = NaN; idx_valid = idx_valid(~isnan(idx_valid));
        gOSI = gOSI(idx_valid); isOS = isOS(idx_valid, :, factor_k);
        idxOS_2 = find(isOS(:, 2) == 1); gOSI_OS = gOSI(idxOS_2);
        gOSI_Tot_k = [gOSI_Tot_k; gOSI]; gOSI_Tot_OS_k = [gOSI_Tot_OS_k; gOSI_OS];
    end
end
end
gOSI_Tot{i} = gOSI_Tot_k; gOSI_Tot_OS{i} = gOSI_Tot_OS_k;
end
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.9, 1]);
for i = 1: length(Dataset_type)
    subplot(2, 4, i); hold on;
    avg1 = mean(gOSI_Tot{i}); plot(avg1 * ones(1, 2), [0 2000], 'k:');
    histogram(gOSI_Tot{i}, x_edge, 'Normalization', 'count');
    axis square; axis([x1, x2, ytickg{i}(1), ytickg{i}(end)]); set(gca, 'XTick', 0: 0.2: 1, 'YTick', ytickg{i});
    title({Dataset_type_title{i}, ['Mean gOSI = ', num2str(avg1, '%.3f')]}, 'FontWeight', 'normal');
    if i == 1, xlabel('gOSI'); ylabel('Number count'); end
    set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
    %
    subplot(2, 4, i + 4); hold on;
    avg2 = mean(gOSI_Tot_OS{i}); plot(avg2 * ones(1, 2), [0 2000], 'k:');
    histogram(gOSI_Tot_OS{i}, x_edge, 'Normalization', 'count');
    axis square; axis([x1, x2, ytickg{i + 4}(1), ytickg{i + 4}(end)]);
    set(gca, 'XTick', 0: 0.2: 1, 'YTick', ytickg{i + 4});
    title({'(OS only)', ['Mean gOSI = ', num2str(avg2, '%.3f')]}, 'FontWeight', 'normal');
    set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
end
%
pause(2); print(gcf, '-dpng', [dir_save, '/gOSI_Overall.png']);
close;


%% Example of tuning curves
rng(1145141919);
load([dir0, '/Analysis_1_Individual_z/L4_cytosolic_GCaMP_Y18_x0y0',...
    '_z320_s2p_NpSize30_Ana1.mat'], 'Tuning_curve_func');
theta_stim_ext_interp = 0: 7.5: 360;
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
% Get OS tuning curves.
for i = 1: length(Dataset_type)
dFF0_mean_OS_o = []; Tuning_par_fit_OS_o = [];
for dataset_k = (Dataset_idx_edge(i) + 1): Dataset_idx_edge(i + 1)
if Dataset_ToUse_SigCorr(dataset_k) == 1
z_value_list = Dataset_z{dataset_k};
factor_k = neuropilFactor_best_idx(dataset_k);
for z_k = 1: length(z_value_list)
    savename = [Dataset_type{i}, '_', Dataset_name{dataset_k}, '_z', num2str(z_value_list(z_k))];
    load([dir0, '/Analysis_1_Individual_z/', savename, '_s2p_NpSize30_Ana1.mat'],...
        'N_neuron', 'dFF0_mean', 'isOS', 'Tuning_par_fit');
    % remove "invalid" neurons
    idx_remove_k = idx_remove{dataset_k}{z_k};    % Comes from neuropil factor auto. No difference...
    idx_valid = 1: N_neuron; idx_valid(idx_remove_k) = NaN; idx_valid = idx_valid(~isnan(idx_valid));
    dFF0_mean = dFF0_mean(idx_valid, :, factor_k);
    isOS = isOS(idx_valid, :, factor_k);
    Tuning_par_fit = Tuning_par_fit(idx_valid, :, factor_k);
    % pick up OS
    idxOS_2 = find(isOS(:, 2) == 1); dFF0_mean_OS = dFF0_mean(idxOS_2, :);
    Tuning_par_fit_OS = Tuning_par_fit(idxOS_2, :);
    % recording
    dFF0_mean_OS_o = [dFF0_mean_OS_o; dFF0_mean_OS];
    Tuning_par_fit_OS_o = [Tuning_par_fit_OS_o; Tuning_par_fit_OS];
end
end
end
%
N_tot = size(dFF0_mean_OS_o, 1);
%
%close; figure; hold on;
%k = randi(N_tot);
%par = Tuning_par_fit_OS_o(k, 1: 6);
%plot(theta_stim, dFF0_mean_OS_o(k, :));
%plot(theta_stim_ext_interp, Tuning_curve_func(theta_stim_ext_interp * pi/180,...
%    par(1), par(2), par(3), par(4), par(5), par(6)));
%plot(theta_stim_ext_interp, interp1(theta_stim, dFF0_mean_OS_o(k, :), theta_stim_ext_interp, 'linear'), 'm:');
%
% We shift the raw shitty tuning curve (but linear interpolated so it can be translated)
i_ctr = find(abs(theta_stim_ext_interp - 180) < eps);
y_tot =  zeros(N_tot, length(theta_stim_ext_interp));
y_avg = zeros(1, length(theta_stim_ext_interp));
for k = 1: N_tot
    y = dFF0_mean_OS_o(k, :);
    y = interp1([theta_stim 360], [y y(1)], theta_stim_ext_interp, 'linear');
    [~, j] = max(y); y = circshift(y, [1, i_ctr - j]);
    y = y / max(y);
    y_tot(k, :) = y; y_avg = y_avg + y / N_tot;
end
% theta_pref, sigma_pref, sigma_oppo, R_pref, R_oppo, R_offset
Err = @(par) Tuning_curve_func((theta_stim_ext_interp(1: 4: end) - 180) * pi/180,...
    par(1), par(2), par(3), par(4), par(5), par(6)) - y_avg(1: 4: end);
IniVal = [0, pi/4, pi/4, 1, 0.5, 0]; par_lb = [- 2 * pi, 0, 0, -100, -100, -100]; par_ub = [2 * pi, pi, pi, 100, 100, 100];
options = optimoptions('lsqnonlin', 'Display', 'none', 'MaxFunEvals', 1200, 'MaxIter', 1200);
par = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
y_avg_fit = Tuning_curve_func((theta_stim_ext_interp - 180) * pi/180,...
    par(1), par(2), par(3), par(4), par(5), par(6));
%
subplot(2, 4, i); hold on;
N_plot = 30; idx = randperm(N_tot); idx = idx(1: N_plot);
for k = 1: N_plot
    plot(theta_stim_ext_interp(1: 4: end) - 180, y_tot(idx(k), 1: 4: end), 'Marker', '.', 'Color', 'b', 'LineWidth', 0.25);
end
%plot(theta_stim_ext_interp - 180, y_avg, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3);
plot(theta_stim_ext_interp - 180, y_avg_fit, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3);
axis([-180 180 0 1.01]); set(gca, 'XTick', -180: 45: 180, 'YTick', 0: 0.25: 1); axis square;
if i == 1, xlabel('\Delta pref. direction'); ylabel('Normalized \DeltaF/F_0'); end
title([Dataset_type_title{i}, ' (OS only)'], 'FontWeight', 'normal');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
subplot(2, 4, i + 4); hold on;
for k = 1: N_plot
    par = Tuning_par_fit_OS_o(idx(k), 1: 6);
    y = Tuning_curve_func((theta_stim_ext_interp - 180) * pi/180,...
        0, par(2), par(3), par(4), par(5), par(6)); y = y / max(y);
    plot(theta_stim_ext_interp - 180, y, 'Marker', '.', 'Color', 'b', 'LineWidth', 0.25);
end
plot(theta_stim_ext_interp - 180, y_avg_fit, 'Color', 'r', 'LineStyle', '--', 'LineWidth', 3);
axis([-180 180 0 1.01]); set(gca, 'XTick', -180: 45: 180, 'YTick', 0: 0.25: 1); axis square;
if i == 1, title('(Gaussian fitting)', 'FontWeight', 'normal'); end
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
end
%
pause(2); print(gcf, '-dpng', [dir_save, '/Tuning_curves_Overall.png']);
close;


