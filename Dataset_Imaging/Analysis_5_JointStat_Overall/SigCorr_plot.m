if ispc, dirD = 'D:'; elseif isunix, dirD = '/media/DATA1'; end
addpath(genpath([dirD, '/Study/CompNeuro/Projects/Functions_simul/']));
dir_data = [dirD, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging'];
load([dir_data, '/General_information.mat']);

d_BinStep_overlap = 5;
%d_BinStep_overlap = 2.5;
d_start = 7.5;
%
Dataset_type_short = {'L2/3 cyt.', 'L2/3 nuc.', 'L4 cyt.', 'L4 soma'};
dir_save = [dir_data, '/Analysis_5_JointStat_Overall/d_Bin', num2str(d_BinStep_overlap)];
%
load([dir_data, '/Analysis_3_JointStat_EachDataset/Results_Datasets_dBin_7.5_',...
    num2str(d_BinStep_overlap, '%.1f'), '.mat']);
load([dir_data, '/Analysis_5_JointStat_Overall/Results_Overall_dBin_7.5_',...
    num2str(d_BinStep_overlap, '%.1f'), '.mat']);
%
idx_start = find(abs(d_BinCenter_overlap - d_start) < 1e-3);
%
Exp1 = @(x, A, lambda, b) A * exp(- x / lambda) + b;
x1 = linspace(5, 50, 101);
Exp2 = @(x, A, lambda, b) A * exp(- (x / lambda).^2) + b;
x2 = linspace(0, 50, 101);
IniVal = [1, 20, 0]; par_lb = [0, 0, 0]; par_ub = [1000, 1000, 1000];
options = optimoptions('lsqnonlin', 'Display', 'none', 'MaxFunEvals', 1200, 'MaxIter', 1200);


d_BinCenter_data4 = cell(1, length(Dataset_type));
SigCorr_mean_data4 = cell(1, length(Dataset_type));
SigCorr_se_data4 = cell(1, length(Dataset_type));
fitpar_exp1_data4 = zeros(length(Dataset_type), 3); fitpar_exp1_err_data4 = zeros(length(Dataset_type), 3);
fitpar_exp2_data4 = zeros(length(Dataset_type), 3); fitpar_exp2_err_data4 = zeros(length(Dataset_type), 3);
%
for i = 1: length(Dataset_type)
    if (i == 1) | (i == 4)
        [~, idx_end] = min(abs(d_BinCenter_overlap - 500));
    else
        [~, idx_end] = min(abs(d_BinCenter_overlap - 240));
    end
    %
    x_data_plot = d_BinCenter_overlap(idx_start: idx_end);
    y_data_plot = SigCorr_mean_overall(idx_start: idx_end, i);
    y_data_se_plot = SigCorr_se_overall(idx_start: idx_end, i);
    idx_nonempty = find(y_data_se_plot > eps);    % = 0, i.e. empty. This is possible after idx_start.
    if length(idx_nonempty) ~= length(x_data_plot), x_data_plot = x_data_plot(idx_nonempty);
    y_data_plot = y_data_plot(idx_nonempty); y_data_se_plot = y_data_se_plot(idx_nonempty); end
    %
    d_BinCenter_data4{i} = x_data_plot; SigCorr_mean_data4{i} = y_data_plot; SigCorr_se_data4{i} = y_data_se_plot;
    %
    x_data_fitting = d_BinCenter_overlap(idx_start: idx_end);
    y_data_fitting = SigCorr_mean_overall(idx_start: idx_end, i);
    y_data_se_fitting = SigCorr_se_overall(idx_start: idx_end, i);
    idx_nonempty = find(y_data_se_fitting > eps);    % = 0, i.e. empty. This is possible after idx_start.
    if length(idx_nonempty) ~= length(x_data_fitting), x_data_fitting = x_data_fitting(idx_nonempty);
    y_data_fitting = y_data_fitting(idx_nonempty); y_data_se_fitting = y_data_se_fitting(idx_nonempty); end
    %
    Err1 = @(par) (Exp1(x_data_fitting, par(1), par(2), par(3)) - y_data_fitting) ./ y_data_se_fitting;
    [fitpar_exp1_data4(i, :), ~, residual, ~, ~, ~, Jacobian] = lsqnonlin(Err1, IniVal, par_lb, par_ub, options);
    CI1 = nlparci(fitpar_exp1_data4(i, :), residual, 'jacobian', Jacobian);
    fitpar_exp1_err_data4(i, :) = ((CI1(:, 2) - CI1(:, 1)) / 2)';
    %
    Err2 = @(par) (Exp2(x_data_fitting, par(1), par(2), par(3)) - y_data_fitting) ./ y_data_se_fitting;
    [fitpar_exp2_data4(i, :), ~, residual, ~, ~, ~, Jacobian] = lsqnonlin(Err2, IniVal, par_lb, par_ub, options);
    CI2 = nlparci(fitpar_exp1_data4(i, :), residual, 'jacobian', Jacobian);
    fitpar_exp2_err_data4(i, :) = ((CI2(:, 2) - CI2(:, 1)) / 2)';
end
Clb_4 = [1 0.5 0; 1 0 0; 0 0 1; 0.5 0 1];



d_BinCenter_data2 = [];
SigCorr_mean_data2 = cell(1, 2); SigCorr_se_data2 = cell(1, 2);
fitpar_exp1_data2 = zeros(2, 3); fitpar_exp1_err_data2 = zeros(2, 3);
fitpar_exp2_data2 = zeros(2, 3); fitpar_exp2_err_data2 = zeros(2, 3);
%
[~, idx_end] = min(abs(d_BinCenter_overlap - 240));
x_data_plot = d_BinCenter_overlap(idx_start: idx_end);
d_BinCenter_data2 = x_data_plot;
x_data_fitting = d_BinCenter_overlap(idx_start: idx_end);
for i = 1: 2
    y_data_plot = SigCorr_mean_layer(idx_start: idx_end, i);
    y_data_se_plot = SigCorr_se_layer(idx_start: idx_end, i);
    idx_nonempty = find(y_data_se_plot > eps);    % = 0, i.e. empty. This is possible after idx_start.
    if length(idx_nonempty) ~= length(x_data_plot), x_data_plot = x_data_plot(idx_nonempty);
    y_data_plot = y_data_plot(idx_nonempty); y_data_se_plot = y_data_se_plot(idx_nonempty); end
    %
    SigCorr_mean_data2{i} = y_data_plot; SigCorr_se_data2{i} = y_data_se_plot;
    %
    y_data_fitting = SigCorr_mean_layer(idx_start: idx_end, i);
    y_data_se_fitting = SigCorr_se_layer(idx_start: idx_end, i);
    idx_nonempty = find(y_data_se_fitting > eps);    % = 0, i.e. empty. This is possible after idx_start.
    if length(idx_nonempty) ~= length(x_data_fitting), x_data_fitting = x_data_fitting(idx_nonempty);
    y_data_fitting = y_data_fitting(idx_nonempty); y_data_se_fitting = y_data_se_fitting(idx_nonempty); end
    %
    Err1 = @(par) (Exp1(x_data_fitting, par(1), par(2), par(3)) - y_data_fitting) ./ y_data_se_fitting;
    [fitpar_exp1_data2(i, :), ~, residual, ~, ~, ~, Jacobian] = lsqnonlin(Err1, IniVal, par_lb, par_ub, options);
    CI1 = nlparci(fitpar_exp1_data4(i, :), residual, 'jacobian', Jacobian);
    fitpar_exp1_err_data2(i, :) = ((CI1(:, 2) - CI1(:, 1)) / 2)';
    %
    Err2 = @(par) (Exp2(x_data_fitting, par(1), par(2), par(3)) - y_data_fitting) ./ y_data_se_fitting;
    [fitpar_exp2_data2(i, :), ~, residual, ~, ~, ~, Jacobian] = lsqnonlin(Err2, IniVal, par_lb, par_ub, options);
    CI2 = nlparci(fitpar_exp1_data4(i, :), residual, 'jacobian', Jacobian);
    fitpar_exp2_err_data2(i, :) = ((CI2(:, 2) - CI2(:, 1)) / 2)';
end
Clb_2 = [1 0 0; 0 0 1];

save([dir_data, '/Analysis_5_JointStat_Overall/SigCorr_plot_d', num2str(d_BinStep_overlap), '.mat'],...
    'd_BinStep_overlap', 'd_BinCenter_data4', 'SigCorr_mean_data4', 'SigCorr_se_data4',...
    'd_BinCenter_data2', 'SigCorr_mean_data2', 'SigCorr_se_data2',...
    'x1', 'Exp1', 'x2', 'Exp2', 'Clb_4', 'Clb_2',...
    'fitpar_exp1_data4', 'fitpar_exp1_err_data4', 'fitpar_exp2_data4', 'fitpar_exp2_err_data4',...
    'fitpar_exp1_data2', 'fitpar_exp1_err_data2', 'fitpar_exp2_data2', 'fitpar_exp2_err_data2');


% % Example of how to use
% Exp_type = 2;
% 
% if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
% load([dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging/',...
%     'Analysis_5_JointStat_Overall/SigCorr_plot_d2.5.mat'],...
%     'd_BinStep_overlap', 'd_BinCenter_data4', 'SigCorr_mean_data4', 'SigCorr_se_data4',...
%     'd_BinCenter_data2', 'SigCorr_mean_data2', 'SigCorr_se_data2',...
%     'x1', 'Exp1', 'x2', 'Exp2', 'Clb_4', 'Clb_2',...
%     'fitpar_exp1_data4', 'fitpar_exp1_err_data4', 'fitpar_exp2_data4', 'fitpar_exp2_err_data4',...
%     'fitpar_exp1_data2', 'fitpar_exp1_err_data2', 'fitpar_exp2_data2', 'fitpar_exp2_err_data2');
% 
% Data_type_title = {'L2/3\_cytosolic\_GCaMP', 'L2/3\_nuclear\_GCaMP', 'L4\_cytosolic\_GCaMP', 'L4\_soma\_GCaMP'};
% Layer_type_title = {'Layer 2/3 (average of all)', 'Layer 4 (average of all)'};
% %
% lwdth = 1.5; mksize = 12.5; cpsize = 7.5; txtsz = 15;
% xylim = [0 50 -0.1 1.01]; xtickc = [0: 5: 50]; ytickc = -0.1: 0.1: 1;
% %
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
% %
% subplot(1, 2, 1); hold on;
% title({'Signal Correlation ~ Horizontal Cortical distance',...
%         ['(Bin Size = ', num2str(d_BinStep_overlap), ' \mum)']}, 'FontWeight', 'normal');
% l = zeros(1, 8); lgdtxt = cell(1, 8);
% for k = 1: 4
%     l(2 * k - 1) = errorbar(d_BinCenter_data4{k}, SigCorr_mean_data4{k}, SigCorr_se_data4{k},...
%         'Color', Clb_4(k, :), 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
%     lgdtxt{2 * k - 1} = Data_type_title{k};
%     if Exp_type == 1, par = fitpar_exp1_data4; par_err = fitpar_exp1_err_data4; 
%         l(2 * k) = plot(x1, Exp1(x1, par(k, 1), par(k, 2), par(k, 3)), 'LineStyle', '-.', 'Color', Clb_4(k, :));
%     elseif Exp_type == 2, par = fitpar_exp2_data4; par_err = fitpar_exp2_err_data4; 
%         l(2 * k) = plot(x2, Exp2(x2, par(k, 1), par(k, 2), par(k, 3)), 'LineStyle', '-.', 'Color', Clb_4(k, :));
%     end
%     lgdtxt{2 * k} = ['(\lambda = ', num2str(par(k, 2), '%.2f'), ' \pm ', num2str(par_err(k, 2), '%.2f'),...
%         ' \mum , A = ', num2str(par(k, 1), '%.3f'), ', b = ', num2str(par(k, 3), '%.3f'), ')'];
% end
% plot([0, xylim(2)], [0 0], 'k--');
% legend(l, lgdtxt, 'FontSize', txtsz);
% %
% subplot(1, 2, 2); hold on;
% l = zeros(1, 2); lgdtxt = cell(1, 2);
% for k = 1: 2
%     l(2 * k - 1) = errorbar(d_BinCenter_data2, SigCorr_mean_data2{k}, SigCorr_se_data2{k},...
%         'Color', Clb_2(k, :), 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
%     lgdtxt{2 * k - 1} = Layer_type_title{k};
%     if Exp_type == 1, par = fitpar_exp1_data2; par_err = fitpar_exp1_err_data2; 
%         l(2 * k) = plot(x1, Exp1(x1, par(k, 1), par(k, 2), par(k, 3)), 'LineStyle', '-.', 'Color', Clb_2(k, :));
%     elseif Exp_type == 2, par = fitpar_exp2_data2; par_err = fitpar_exp2_err_data2; 
%         l(2 * k) = plot(x2, Exp2(x2, par(k, 1), par(k, 2), par(k, 3)), 'LineStyle', '-.', 'Color', Clb_2(k, :));
%     end
%     lgdtxt{2 * k} = ['(\lambda = ', num2str(par(k, 2), '%.2f'), ' \pm ', num2str(par_err(k, 2), '%.2f'),...
%         ' \mum , A = ', num2str(par(k, 1), '%.3f'), ', b = ', num2str(par(k, 3), '%.3f'), ')'];
% end
% plot([0, xylim(2)], [0 0], 'k--');
% legend(l, lgdtxt, 'FontSize', txtsz);
% %
% for i = 1: 2
%     subplot(1, 2, i);
%     axis(xylim); set(gca, 'XTick', xtickc, 'YTick', ytickc);
%     axis square; grid on; 
%     xlabel('Horizontal Cortical Distance (\mum)');
%     set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
% end

