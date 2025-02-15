if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
addpath(genpath([dir_0, '/Study/CompNeuro/Projects/Functions_simul/']));
cd([dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Model_SigCorr/LR_continuum/Core']);

load([dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging/',...
    'Analysis_5_JointStat_Overall/SigCorr_plot_d5.mat'],...
    'fitpar_exp2_data2', 'SigCorr_mean_data2', 'SigCorr_se_data2', 'Exp2', 'Exp1');
fitpar_exp2_L4 = fitpar_exp2_data2(2, :);    % L4 input, use exp2 model
%
d_micron = 0: 0.1: 750;
W = [3.2, 1.08, -3.5; 6.4, 0.8, -6.15];
sigma_micron = [144, 146, 110; 3.75, 20, 9.25];
kappa = [0.075, 0.125, 0.125];
%
[FR_n_L23E, FR_n_L23I, Cov_r_L23E, Cov_r_L23I] =...
    LR_continuum_FT_func(d_micron, fitpar_exp2_L4, W, sigma_micron, kappa);
%
options = optimoptions('lsqnonlin', 'Display', 'none', 'MaxFunEvals', 1200,...
    'MaxIter', 1200, 'TolX', 1e-8, 'TolFun', 1e-8);
[~, j] = find(d_micron - [7.5: 5: 52.5]' == 0);
ErrA = @(Amp) (Amp * Cov_r_L23E(j)' - SigCorr_mean_data2{1}(1: 10)) ./ SigCorr_se_data2{1}(1: 10);
A = lsqnonlin(ErrA, 2500, 0, 1e4, options);
%A = 30;
Corr_r_L23E = Cov_r_L23E * A;
%
fitting_exp_type = 2;
IniVal = [Corr_r_L23E(1), 20, 0]; par_lb = [0, 0, 0]; par_ub = [1000, 1000, 1000];
if fitting_exp_type == 1
    Err = @(par) 1e5 * (Exp1(d_micron, par(1), par(2), par(3)) - Corr_r_L23E);
elseif fitting_exp_type == 2
    Err = @(par) 1e5 * (Exp2(d_micron, par(1), par(2), par(3)) - Corr_r_L23E);
end
[fitpar_Corr_r_L23E, ~, residual, ~, ~, ~, Jacobian] = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
CI = nlparci(fitpar_Corr_r_L23E, residual, 'jacobian', Jacobian);
fitpar_err_Corr_r_L23E = ((CI(:, 2) - CI(:, 1)) / 2)';




load([dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging/',...
    'Analysis_5_JointStat_Overall/SigCorr_plot_d5.mat'],...
    'd_BinCenter_data2', 'SigCorr_mean_data2', 'SigCorr_se_data2',...
    'fitpar_exp1_data2', 'fitpar_exp1_err_data2', 'Exp1',...
    'fitpar_exp2_data2', 'fitpar_exp2_err_data2', 'Exp2', 'Clb_2');
x1 = 7.5: 0.1: 50; x2 = 0: 0.1: 50;
%
% plot order: L4 data, L2/3 data, L2/3 model
N_legend = 3; legend_order_data = [1 2]; legend_order_model = [3];
l = zeros(1, N_legend); lgdtxt = cell(1, N_legend);
%
% Data
data_type_title = {'L2/3 data', 'L4 data'};
lwdth = 1.5; mksize = 12.5; cpsize = 7.5; txtsz = 15;
xylim = [0 50 -0.1 1.01]; xtickc = [0: 5: 50]; ytickc = -0.1: 0.1: 1;
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.55, 1]); hold on;
for k = 1: length(data_type_title)
    m = legend_order_data(k);
    l(m) = errorbar(d_BinCenter_data2, SigCorr_mean_data2{k}, SigCorr_se_data2{k},...
        'Color', Clb_2(k, :), 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
    if fitting_exp_type == 1
        par = fitpar_exp1_data2; par_err = fitpar_exp1_err_data2;
        plot(x1, Exp1(x1, par(k, 1), par(k, 2), par(k, 3)),...
            'LineStyle', '--', 'LineWidth', lwdth, 'Color', Clb_2(k, :));
    elseif fitting_exp_type == 2
        par = fitpar_exp2_data2; par_err = fitpar_exp2_err_data2;
        plot(x2, Exp2(x2, par(k, 1), par(k, 2), par(k, 3)),...
            'LineStyle', '--', 'LineWidth', lwdth, 'Color', Clb_2(k, :));
    end
    lgdtxt{m} = [data_type_title{k}, ' (A = ', num2str(par(k, 1), '%.2f'),...
        ', \lambda = ', num2str(par(k, 2), '%.2f'), ' \pm ', num2str(par_err(k, 2), '%.2f'), ' \mum)'];
end
plot([0, xylim(2)], [0 0], 'k--');
axis(xylim); set(gca, 'XTick', xtickc, 'YTick', ytickc);
axis square; grid on;
xlabel('Horizontal Cortical Distance (\mum)');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
% Model
idx_bound = find(d_micron == 7.5); Clb_model = [0.5, 0, 0];
plot(d_micron(1: idx_bound - 1), Corr_r_L23E(1: idx_bound - 1),...
    'Color', Clb_model, 'LineStyle', '--', 'LineWidth', lwdth);
l(legend_order_model) = plot(d_micron(idx_bound: end), Corr_r_L23E(idx_bound: end),...
    'Color', Clb_model, 'LineWidth', lwdth + 1);
lgdtxt{legend_order_model} =...
    ['L2/3 model (continuum) (A = ', num2str(fitpar_Corr_r_L23E(1), '%.2f'),...
        ', \lambda = ', num2str(fitpar_Corr_r_L23E(2), '%.2f'), ' \pm ',...
        num2str(fitpar_err_Corr_r_L23E(2), '%.2f'), ' \mum)'];
%
legend(l, lgdtxt);
text_sigmaB = ['\sigmaB = [', num2str(sigma_micron(1, 1)), ', ',...
    num2str(sigma_micron(1, 2)), ', ', num2str(sigma_micron(1, 3)), ']'];
text_sigmaN = ['\sigmaN = [', num2str(sigma_micron(2, 1)), ', ',...
    num2str(sigma_micron(2, 2)), ', ', num2str(sigma_micron(2, 3)), ']'];
text_kappa = ['\kappa = [', num2str(kappa(1)), ', ', num2str(kappa(2)), ', ', num2str(kappa(3)), ']'];
text_W = ['W = [', num2str(W(1, 1)), ', ', num2str(W(1, 2)), ', ', num2str(W(1, 3)),...
    '; ', num2str(W(2, 1)), ', ', num2str(W(2, 2)), ', ', num2str(W(2, 3)), ']'];
suptitle({[text_sigmaB, '; ', text_sigmaN, ' (\mum)'], [text_kappa, '; ', text_W, ';']}, 6, 0.95);
%
pause(2); print(gcf, '-dpng', 'Test_LR_continuum.png'); close;

