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
%
if d_BinStep_overlap == 5
    idx_start = find(abs(d_BinCenter_overlap - 7.5) < 1e-3);
elseif d_BinStep_overlap == 2.5
    idx_start = find(abs(d_BinCenter_overlap - 7.5) < 1e-3);
end
%
Dataset_type_title = {'L2/3\_cytosolic\_GCaMP', 'L2/3\_nuclear\_GCaMP', 'L4\_cytosolic\_GCaMP', 'L4\_soma\_GCaMP'};
Dataset_name_title = {'KF19', 'KF20', 'D12', 'D3', 'D4', 'Y18\_x0y0', 'Y18\_x-256y55',...
    'Y22\_x0y0', 'Y22\_x-401y418', 'Y24\_x0y0', 'Y24\_x426y133', 'Y25\_x0y0', 'Y25\_x10y-256', 'Y26', 'M199'};
%
Exp1 = @(x, A, lambda, b) - A * exp(- x / lambda) + b;
x1 = linspace(2.5, 50, 101);
Exp2 = @(x, A, lambda, b) - A * exp(- (x / lambda).^2) + b;
x2 = linspace(0, 50, 101);
IniVal = [45, 10, 45]; par_lb = [0, 1e-5, 0]; par_ub = [1000, 1000, 100];
options = optimoptions('lsqnonlin', 'Display', 'none', 'MaxFunEvals', 1200, 'MaxIter', 1200);
%
lwdth = 1.5; mksize = 12.5; cpsize = 7.5; txtsz = 15;
Clb_4 = [1 0.5 0; 1 0 0; 0 0 1; 0.5 0 1];
xylim = [0 50 0 55]; xtickc = [0: 5: 50]; ytickc = 0: 5: 55;


%% Overlap results of each dataset
%
for i = 1: length(Dataset_type)
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.55, 1]); hold on;
suptitle({['\Delta\theta\_pref. ~ Horizontal Cortical distance, ', Dataset_type_title{i}],...
        ['(Bin Size = ', num2str(d_BinStep_overlap), ' \mum)']}, 8, 0.95);
l = zeros(1, Dataset_N(i)); lgdtxt = cell(1, Dataset_N(i));
clb = jet(Dataset_N(i)) * 0.875; clb = clb(end: -1: 1, :);
for k = 1: Dataset_N(i)
    dataset_k = Dataset_idx_edge(i) + k;
    factor_k = neuropilFactor_best_idx(dataset_k);
    %
    if ROI_size_Tot(dataset_k, 1) < 300    % 256
        [~, idx_end] = min(abs(d_BinCenter_overlap - 240));
    else
        [~, idx_end] = min(abs(d_BinCenter_overlap - 500));
    end
    %
    x_data = d_BinCenter_overlap(idx_start: idx_end);
    y_data = dtheta_pref_mean_datasets(idx_start: idx_end, dataset_k, factor_k);
    y_data_se = dtheta_pref_se_datasets(idx_start: idx_end, dataset_k, factor_k);
    idx_nonempty = find(y_data_se > eps);    % = 0, i.e. empty. This is possible after idx_start.
    if length(idx_nonempty) ~= length(x_data)
        x_data = x_data(idx_nonempty); y_data = y_data(idx_nonempty); y_data_se = y_data_se(idx_nonempty);
    end
    if y_data(2) < y_data(1), x_data = x_data(2: end); y_data = y_data(2: end); y_data_se = y_data_se(2: end); end
    if (Dataset_ToUse_dtheta(dataset_k) == 0), lst = '--'; else, lst = '-'; end
    l(k) = errorbar(x_data, y_data, y_data_se, 'Color', clb(k, :), 'LineStyle', lst,...
        'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
    lgdtxt{k} = Dataset_name_title{dataset_k};
end
plot([0, xylim(2)], [45 45], 'k--');
axis(xylim); set(gca, 'XTick', xtickc, 'YTick', ytickc);
axis square; grid on; legend(l, lgdtxt, 'FontSize', txtsz, 'Location', 'southeast');
xlabel('Horizontal Cortical Distance (\mum)');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
pause(2); print(gcf, '-dpng', [dir_save, '/dtheta_pref_', Dataset_type{i}, '.png']);
close;
end


%% Overall results of L2/3c, L2/3n, L4c
for k = 1: 2
fitpar = zeros(3, 3); fitpar_err = zeros(3, 3);    % row: par; column: lb / ub.
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.55, 1]); hold on;
suptitle({'\Delta\theta\_pref. ~ Horizontal Cortical distance',...
        ['(Bin Size = ', num2str(d_BinStep_overlap), ' \mum)']}, 8, 0.95);
l = zeros(1, 2 * length(Dataset_N)); lgdtxt = cell(1, 2 * length(Dataset_N));
for i = 1: length(Dataset_type)
    if (i == 1) | (i == 4)
        [~, idx_end] = min(abs(d_BinCenter_overlap - 500));
    else
        [~, idx_end] = min(abs(d_BinCenter_overlap - 240));
    end
    %
    x_data = d_BinCenter_overlap(idx_start: idx_end);
    y_data = dtheta_pref_mean_overall(idx_start: idx_end, i);
    y_data_se = dtheta_pref_se_overall(idx_start: idx_end, i);
    idx_nonempty = find(y_data_se > eps);    % = 0, i.e. empty. This is possible after idx_start.
    if length(idx_nonempty) ~= length(x_data)
        x_data = x_data(idx_nonempty); y_data = y_data(idx_nonempty); y_data_se = y_data_se(idx_nonempty);
    end
    if y_data(2) < y_data(1), x_data = x_data(2: end); y_data = y_data(2: end); y_data_se = y_data_se(2: end); end
    l(2 * (i - 1) + 1) = errorbar(x_data, y_data, y_data_se, 'Color', Clb_4(i, :),...
        'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
    lgdtxt{2 * (i - 1) + 1} = Dataset_type_title{i};
    %
    if k == 1
        Err = @(par) (Exp1(x_data, par(1), par(2), par(3)) - y_data) ./ y_data_se;
    elseif k == 2
        Err = @(par) (Exp2(x_data, par(1), par(2), par(3)) - y_data) ./ y_data_se;
    end
    [fitpar(i, :), ~, residual, exitflag, ~, ~, Jacobian] = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
    CI = nlparci(fitpar(i, :), residual, 'jacobian', Jacobian); fitpar_err(i, :) = ((CI(:, 2) - CI(:, 1)) / 2)';
    if k == 1
        l(2 * i) = plot(x1, Exp1(x1, fitpar(i, 1), fitpar(i, 2), fitpar(i, 3)), 'LineStyle', '-.', 'Color', Clb_4(i, :));
    elseif k == 2
        l(2 * i) = plot(x2, Exp2(x2, fitpar(i, 1), fitpar(i, 2), fitpar(i, 3)), 'LineStyle', '-.', 'Color', Clb_4(i, :));
    end
    lgdtxt{2 * i} = ['(\lambda = ', num2str(fitpar(i, 2), '%.2f'), ' \pm ', num2str(fitpar_err(i, 2), '%.2f'),...
        ' \mum , A = ', num2str(fitpar(i, 1), '%.3f'), ', b = ', num2str(fitpar(i, 3), '%.3f'), ')'];
end
plot([0, xylim(2)], [45 45], 'k--');
axis(xylim); set(gca, 'XTick', xtickc, 'YTick', ytickc);
axis square; grid on; legend(l, lgdtxt, 'FontSize', txtsz, 'Location', 'southeast');
xlabel('Horizontal Cortical Distance (\mum)');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
pause(2); print(gcf, '-dpng', [dir_save, '/dtheta_pref_Overall_Exp', num2str(k), '.png']);
close;
end

end

