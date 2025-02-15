if ispc, dirD = 'D:'; elseif isunix, dirD = '/media/DATA1'; end
addpath(genpath([dirD, '/Study/CompNeuro/Projects/Functions_simul/']));
dir0 = [dirD, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging'];
load([dir0, '/General_information.mat']);

dir_name = [dir0, '/Analysis_4_NeuropilFactor_EachDataset/d_Bin5'];
if ~exist(dir_name, 'dir'), mkdir(dir_name); end
dir_name = [dir0, '/Analysis_4_NeuropilFactor_EachDataset/d_Bin2.5'];
if ~exist(dir_name, 'dir'), mkdir(dir_name); end
dir_name = [dir0, '/Analysis_4_NeuropilFactor_EachDataset/d_Bin5/PrefOriMap'];
if ~exist(dir_name, 'dir'), mkdir(dir_name); end
dir_name = [dir0, '/Analysis_4_NeuropilFactor_EachDataset/d_Bin5/PrefOriDistribution'];
if ~exist(dir_name, 'dir'), mkdir(dir_name); end
dir_name = [dir0, '/Analysis_4_NeuropilFactor_EachDataset/d_Bin5/SigCorr'];
if ~exist(dir_name, 'dir'), mkdir(dir_name); end
dir_name = [dir0, '/Analysis_4_NeuropilFactor_EachDataset/d_Bin2.5/PrefOriMap'];
if ~exist(dir_name, 'dir'), mkdir(dir_name); end
dir_name = [dir0, '/Analysis_4_NeuropilFactor_EachDataset/d_Bin2.5/PrefOriDistribution'];
if ~exist(dir_name, 'dir'), mkdir(dir_name); end
dir_name = [dir0, '/Analysis_4_NeuropilFactor_EachDataset/d_Bin2.5/SigCorr'];
if ~exist(dir_name, 'dir'), mkdir(dir_name); end


Dataset_type_title = {'L2/3\_cytosolic\_GCaMP', 'L2/3\_nuclear\_GCaMP', 'L4\_cytosolic\_GCaMP',  'L4\_soma\_GCaMP'};
Dataset_name_title = {'KF19', 'KF20', 'D12', 'D3', 'D4', 'Y18\_x0y0', 'Y18\_x-256y55',...
    'Y22\_x0y0', 'Y22\_x-401y418', 'Y24\_x0y0', 'Y24\_x426y133', 'Y25\_x0y0', 'Y25\_x10y-256', 'Y26', 'M199'};
%
d_BinStep_individual = 7.5;
for d_BinStep_overlap = [2.5 5]
load([dir0, '/Analysis_3_JointStat_EachDataset/Results_Datasets_dBin_',...
    num2str(d_BinStep_individual, '%.1f'), '_', num2str(d_BinStep_overlap, '%.1f'), '.mat']);
neuron_idx_dataset = [zeros(1, N_neuropilFactor); cumsum(N_neuron_OS_2_datasets, 1)];
if d_BinStep_overlap == 5
    idx_start = find(abs(d_BinCenter_overlap - 7.5) < 1e-3);
elseif d_BinStep_overlap == 2.5
    idx_start = find(abs(d_BinCenter_overlap - 7.5) < 1e-3);
end
%
plot_idx = [1, 2, 3, 5, 6, 7, 8];
hsv_clb = hsv * 0.875; mksize = 10; %scatter_edge_width = 0.25;
%
theta1 = 0; theta2 = 180; dtheta = 7.5; theta_edge = (theta1 - dtheta / 2): dtheta: (theta2 + dtheta / 2);
%
Clb = [[0: 0.25: 1; zeros(1, 5); 1: -0.25: 0]'; [1 0.6 0.6]; [0.5 0.5 0.5]];
lwdth = 1.5; mksize = 12.5; cpsize = 7.5; txtsz = 15;
xylim = [0 200 -0.1 1]; xtickc = [0: 10: 200]; ytickc = -0.1: 0.05: 1;
 

xticklc = {'0', '', '', '', '', '50', '', '', '', '', '100',...
    '', '', '', '', '150', '', '', '', '', '200'};
yticklc = {'-0.1', '', '0', '', '0.1', '', '0.2', '', '0.3', '', '0.4', '',...
    '0.5', '', '0.6', '', '0.7', '', '0.8', '', '0.9', '', '1'};
%
Exp12 = @(x, An, ln, Ab, lb, b) An * exp(- x / ln) + Ab * exp(- x / lb) + b;
IniVal = [repmat([1, 10, 0.1, 0.5], 1, 3), [1, 10, 0.1, 0.1], repmat([1, 10, 0.1, 0.01], 1, 2), 100];
par_lb = 0 * ones(1, 25); par_ub = 500 * ones(1, 25);
options = optimoptions('lsqnonlin', 'Display', 'none', 'MaxFunEvals', 1200, 'MaxIter', 1200);




for i = 1: length(Dataset_type)
for dataset_k = (Dataset_idx_edge(i) + 1): Dataset_idx_edge(i + 1)
savename0 = [dir0, '/Analysis_4_NeuropilFactor_EachDataset/d_Bin', num2str(d_BinStep_overlap)];
savename1 = [Dataset_type{i}, '_', Dataset_name{dataset_k}];
%
%% Pref. orientation maps
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 0.85]);
suptitle([Dataset_type_title{i}, ', ', Dataset_name_title{dataset_k}], 4, 0.975);
for factor_k = 1: N_neuropilFactor
idx_k = [neuron_idx_dataset(dataset_k, factor_k) + 1: neuron_idx_dataset(dataset_k + 1, factor_k)];
ROI_centeroid_k = ROI_centeroid_OS_2_Tot{factor_k}(idx_k, :);
Pref_orientation_k = Pref_orientation_OS_Tot{factor_k}(idx_k);
N_neuron_k = N_neuron_OS_2_datasets(dataset_k, factor_k);
%
color_idx = ceil((Pref_orientation_k / 180) * size(hsv_clb, 1)); color_idx(color_idx == 0) = 1;
color_array = zeros(N_neuron_k, 3);
for neuron_k = 1: N_neuron_k, color_array(neuron_k, :) = hsv_clb(color_idx(neuron_k), :); end
%
subplot(2, 4, plot_idx(factor_k));
scatter(ROI_centeroid_k(:, 2), ROI_size_Tot(dataset_k, 2) - ROI_centeroid_k(:, 1),...    % xy
    mksize, color_array, 'fill', 'MarkerEdgeColor', 'none');%, 'LineWidth', scatter_edge_width); 
axis square; axis([0, ROI_size_Tot(dataset_k, 1), 0, ROI_size_Tot(dataset_k, 2)]);
if ROI_size_Tot(dataset_k, 1) < 300    % 256
    set(gca, 'XTick', 0: 32: 256, 'YTick', 0: 32: 256);
elseif ROI_size_Tot(dataset_k, 1) < 1000    % 1024 * 0.821, 840
    set(gca, 'XTick', 0: 100: 800, 'YTick', 0: 100: 800);
else    % 1024
    set(gca, 'XTick', 0: 100: 1000, 'YTick', 0: 100: 1000);
end
if factor_k ~= N_neuropilFactor
    title(['\alpha = ', num2str(neuropilFactor_manual(factor_k)),...
        ' (N = ', num2str(N_neuron_k), ')'], 'FontWeight', 'normal');
else
    title(['Auto \alpha (avg. ', num2str(neuropilFactor_auto_avg_datasets(dataset_k), '%.2f'),...
        ') (N = ', num2str(N_neuron_k), ')'], 'FontWeight', 'normal');
    colormap(hsv_clb); clb = colorbar; clb.Position = [0.92, 0.11, 0.01, 0.342];
    set(clb, 'YTick', 0: 1/4: 1, 'YTickLabel', {'0^o', '45^o', '90^o', '135^o', '180^o'});
end
if factor_k == 1, xlabel('x (\mum)'); ylabel('y (\mum)'); end
set(gca, 'box', 'on', 'TickDir', 'out');
end
%
pause(2); print(gcf, '-dpng', [savename0, '/PrefOriMap/', savename1, '_PrefOriMap.png']);
close;
%
%% Distribution of Preferred Orientations
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 0.85]);
suptitle([Dataset_type_title{i}, ', ', Dataset_name_title{dataset_k}], 4, 0.975);
for factor_k = 1: N_neuropilFactor
idx_k = [neuron_idx_dataset(dataset_k, factor_k) + 1: neuron_idx_dataset(dataset_k + 1, factor_k)];
Pref_orientation_k = Pref_orientation_OS_Tot{factor_k}(idx_k);
N_neuron_k = N_neuron_OS_2_datasets(dataset_k, factor_k);
%
subplot(2, 4, plot_idx(factor_k));
histogram(Pref_orientation_k, theta_edge, 'Normalization', 'count', 'FaceColor', [0 0.5 1]);
[pval, ~] = circ_otest(Pref_orientation_k * 2 * (pi / 180), dtheta * (pi / 180));
if pval >= 0.01, title2 = ['p = ', num2str(pval, '%.3f')];
else, title2 = ['p < 10^{', num2str(ceil(log10(pval))), '}']; end
axis square; xlim([theta1 theta2]); set(gca, 'XTick', 0: 15: 180, 'XTickLabel',...
    {'0', '', '30', '', '60', '', '90', '', '120', '', '150', '', '180'});
if factor_k ~= N_neuropilFactor
    title(['\alpha = ', num2str(neuropilFactor_manual(factor_k)),...
        ', ', title2, ' (N = ', num2str(N_neuron_k), ')'], 'FontWeight', 'normal');
else
    title(['Auto \alpha (avg. ', num2str(neuropilFactor_auto_avg_datasets(dataset_k), '%.2f'),...
        '), ', title2, ' (N = ', num2str(N_neuron_k), ')'], 'FontWeight', 'normal');
end
if factor_k == 1
    xlabel('Pref. orientation'); ylabel('Number count');
    ax = gca; ylim_tmp = ax.YLim; ytick_tmp = ax.YTick;
else
    ylim(ylim_tmp); set(gca, 'YTick', ytick_tmp);
end
set(gca, 'box', 'off', 'TickDir', 'out');
end
%
pause(2); print(gcf, '-dpng', [savename0, '/PrefOriDistribution/', savename1, '_PrefOriDistribution.png']);
close;
%
%% Signal correlation
% Do parameter fitting for only KF20 / Y18-x0y0 -- will be used in the paper.
if ROI_size_Tot(dataset_k, 1) < 300    % 256
    [~, idx_end] = min(abs(d_BinCenter_overlap - 240));
else
    [~, idx_end] = min(abs(d_BinCenter_overlap - 500));
end
x_data_tot = cell(1, N_neuropilFactor); y_data_tot = cell(1, N_neuropilFactor);
l = zeros(1, N_neuropilFactor); lgdtxt = cell(1, N_neuropilFactor);
%
figure;
if any(dataset_k == [2, 6])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
    suptitle({'Signal Correlation ~ Horizontal Cortical distance',...
        ['(Bin Size = ', num2str(d_BinStep_overlap), ' \mum)']}, 8, 0.965);
    subplot(2, 4, [1 2 5 6]);
else
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.7, 1]);
    suptitle({'Signal Correlation ~ Horizontal Cortical distance',...
        ['(Bin Size = ', num2str(d_BinStep_overlap), ' \mum)']}, 8, 0.95);
end
hold on;
for factor_k = 1: N_neuropilFactor
    x_data = d_BinCenter_overlap(idx_start: idx_end); n1 = length(x_data);
    y_data = SigCorr_mean_datasets(idx_start: idx_end, dataset_k, factor_k);
    y_data_se = SigCorr_se_datasets(idx_start: idx_end, dataset_k, factor_k);
    idx_nonempty = find(y_data_se > eps);    % = 0, i.e. empty. This is possible after idx_start.
    if length(idx_nonempty) ~= length(x_data)
        x_data = x_data(idx_nonempty); y_data = y_data(idx_nonempty); y_data_se = y_data_se(idx_nonempty);
    end
    l(factor_k) = errorbar(x_data, y_data, y_data_se, 'Color', Clb(factor_k, :),...
        'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
    x_data_tot{factor_k} = x_data; y_data_tot{factor_k} = y_data;
    if factor_k ~= N_neuropilFactor
        lgdtxt{factor_k} = ['\alpha = ', num2str(neuropilFactor_manual(factor_k))];
    else
        lgdtxt{factor_k} = ['Auto \alpha (avg. ', num2str(neuropilFactor_auto_avg_datasets(dataset_k), '%.2f'), ')'];
    end
end
plot([0, xylim(2)], [0 0], 'k--');
axis(xylim); set(gca, 'XTick', xtickc, 'YTick', ytickc, 'XTickLabel', xticklc, 'YTickLabel', yticklc);
grid on; legend(l, lgdtxt, 'FontSize', txtsz);
xlabel('Horizontal Cortical Distance (\mum)');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
if any(dataset_k == [2, 6])
% (An1, ln1, Ab1, b1, An2, ln2, Ab2, b2, An3, ln3, Ab3, b3,...
    % An4, ln4, Ab4, b4, An5, ln5, Ab5, b5, An6, ln6, Ab6, b6, lb)
Err = @(par) [Exp12(x_data_tot{1}, par(1), par(2), par(3), par(25), par(4)) - y_data_tot{1};...
    Exp12(x_data_tot{2}, par(5), par(6), par(7), par(25), par(8)) - y_data_tot{2};...
    Exp12(x_data_tot{3}, par(9), par(10), par(11), par(25), par(12)) - y_data_tot{3};...
    Exp12(x_data_tot{4}, par(13), par(14), par(15), par(25), par(16)) - y_data_tot{4};...
    Exp12(x_data_tot{5}, par(17), par(18), par(19), par(25), par(20)) - y_data_tot{5};...
    Exp12(x_data_tot{6}, par(21), par(22), par(23), par(25), par(24)) - y_data_tot{6}];
[par, ~, residual, exitflag, ~, ~, Jacobian] = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
CI = nlparci(par, residual, 'jacobian', Jacobian);
A_narrow = par(1: 4: 21); A_narrow_err = (CI(1: 4: 21, 2)' - CI(1: 4: 21, 1)') / 2;
A_broad = par(3: 4: 23); A_broad_err = (CI(3: 4: 23, 2)' - CI(3: 4: 23, 1)') / 2;
lambda_narrow = par(2: 4: 22); lambda_narrow_err = (CI(2: 4: 22, 2)' - CI(2: 4: 22, 1)') / 2;
b0 = par(4: 4: 24); b0_err = (CI(4: 4: 24, 2)' - CI(4: 4: 24, 1)') / 2;
lambda_broad = par(25); lambda_broad_err = (CI(25, 2)' - CI(25, 1)') / 2;
%
subplot(2, 4, 3); hold on;
errorbar(neuropilFactor_manual, lambda_broad * ones(1, 6), lambda_broad_err * ones(1, 6), 'Color', [1, 0, 0],...
    'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
errorbar(neuropilFactor_manual, lambda_narrow, lambda_narrow_err, 'Color', [0, 0.5, 1],...
    'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
if dataset_k == 2
    axis([-0.01 1.26 -0.5 250]); set(gca, 'XTick', 0: 0.25: 1.25, 'YTick', 0: 50: 250);
elseif dataset_k == 6
    axis([-0.01 1.26 -0.5 30]); set(gca, 'XTick', 0: 0.25: 1.25, 'YTick', 0: 5: 30);
end
title('\lambda (broad & narrow)', 'FontWeight', 'normal');
legend({'\lambda\_broad', '\lambda\_narrow'}, 'Location', 'east');
axis square; grid on; set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
subplot(2, 4, 4); hold on;
errorbar(neuropilFactor_manual, A_broad, A_broad_err, 'Color', [0, 0.5, 1],...
    'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
if dataset_k == 2
    axis([-0.01 1.26 -0.001 0.25]); set(gca, 'XTick', 0: 0.25: 1.25, 'YTick', 0: 0.05: 0.25);
elseif dataset_k == 6
    axis([-0.01 1.26 -0.001 0.75]); set(gca, 'XTick', 0: 0.25: 1.25, 'YTick', 0: 0.15: 0.75);
end
title('A\_broad', 'FontWeight', 'normal');
axis square; grid on; set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
subplot(2, 4, 7); hold on;
errorbar(neuropilFactor_manual, A_narrow, A_narrow_err, 'Color', [0, 0.5, 1],...
    'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
if dataset_k == 2
    axis([-0.01 1.26 -0.005 4]); set(gca, 'XTick', 0: 0.25: 1.25, 'YTick', 0: 1: 4);
elseif dataset_k == 6
    axis([-0.01 1.26 -0.005 2]); set(gca, 'XTick', 0: 0.25: 1.25, 'YTick', 0: 0.5: 2);
end
title('A\_narrow', 'FontWeight', 'normal');
axis square; grid on; set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
subplot(2, 4, 8); hold on;
errorbar(neuropilFactor_manual, b0, b0_err, 'Color', [0, 0.5, 1],...
    'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
if dataset_k == 2
    axis([-0.01 1.26 -0.001 0.6]); set(gca, 'XTick', 0: 0.25: 1.25, 'YTick', 0: 0.15: 0.6);
elseif dataset_k == 6
    axis([-0.01 1.26 -0.0001 0.05]); set(gca, 'XTick', 0: 0.25: 1.25, 'YTick', 0: 0.01: 0.05);
end
title('b', 'FontWeight', 'normal');
axis square; grid on; set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
end
%
pause(2); print(gcf, '-dpng', [savename0, '/SigCorr/', savename1, '_SigCorr.png']);
savename_full = [savename0, '/SigCorr/', savename1, '_SigCorr'];
print(gcf, '-dpng', [savename_full, '.png']);
if any(dataset_k == [2, 3, 6, 15])
    print(gcf, '-depsc2', '-r1500', [savename_full, '.eps']);
end
close;

end
end

end

