if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
addpath(genpath([dir_0, '/Study/CompNeuro/Projects/Functions_simul/']));


for Parset_i = [2 3 5]
%Parset_i = 4;

% Load data
load([dir_0, '/Study/CompNeuro/Projects/Micro-clustering',...
    '/Dataset_Imaging/Noise_correlation/Results_Noise_correlation.mat'],...
    'd_BinCenter', 'Corr_avg', 'Corr_se', 'DualExp2', 'FitPar', 'FitPar_err');
%
% DualExp2 = @(x, A1, A2, lambda1, lambda2, b)...
%     A1 * exp(- (x / lambda1).^2) + A2 * exp(- (x / lambda2).^2) + b;
x2 = 0: 0.1: 100;
%
% Load model
dir_local = [dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Model_NoiseCorr'];
load([dir_local, '/Results/Results_Total_', num2str(Parset_i), '.mat'],...
    'd_BinCenter_F', 'd_BinCenter_E', 'dtheta_bound1', 'dtheta_bound2',...
    'NoiseCorr_mean_F', 'NoiseCorr_se_F', 'fitpar_noi_F', 'fitpar_err_noi_F',...
    'NoiseCorr_mean_E', 'NoiseCorr_se_E', 'fitpar_noi_E', 'fitpar_err_noi_E',...
    'NoiseCorr_mean_E_dtheta1', 'NoiseCorr_se_E_dtheta1', 'fitpar_noi_E_dtheta1', 'fitpar_err_noi_E_dtheta1',...
    'NoiseCorr_mean_E_dtheta2', 'NoiseCorr_se_E_dtheta2', 'fitpar_noi_E_dtheta2', 'fitpar_err_noi_E_dtheta2');
%
load([dir_local, '/Parameters/Parameters_Recurrent_1.mat'], 'param');
T_on = param.T_on; clear param
%
lwdth = 1.5; mksize = 12.5; cpsize = 7.5; txtsz = 12;
xylim = [0 50 -0.02 0.501]; xtickc = [0: 5: 50]; ytickc = 0: 0.05: 0.5;



clear l lgdtxt
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 0.7]);
suptitle(['Noise Correlation (evoked rsp. only), Data (NF = 0.75) v.s. Model'], 6, 0.93);
%
% L4
subplot(1, 3, 1); hold on;
l = zeros(1, 2); Clr_L4 = [0 0.5 1];
%
% Model only
l(1) = errorbar(d_BinCenter_F, NoiseCorr_mean_F, NoiseCorr_se_F,...
    'Color', Clr_L4, 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
lgdtxt{1} = 'L4 model (all pairs)';
%
par = fitpar_noi_F; par_err = fitpar_err_noi_F;
l(2) = plot(x2, DualExp2(x2, par(1), 0, par(3), 114514, 0),...
    'LineStyle', '--', 'LineWidth', lwdth, 'Color', Clr_L4);
lgdtxt{2} =  ['(An = ', num2str(par(1), '%.2f'), ', \lambdan = ',...
    num2str(par(3), '%.2f'), ' \pm ', num2str(par_err(3), '%.2f'), ' \mum)'];

%
plot([0, xylim(2)], [0 0], 'k--');
legend(l, lgdtxt, 'FontSize', txtsz);
title({'Layer 4', ['(Time bin ', num2str(T_on), ' ms)']}, 'FontWeight', 'normal');
%
%axis(xylim); set(gca, 'XTick', xtickc, 'YTick', ytickc);
%axis square; grid on;
%xlabel('Horizontal Cortical Distance (\mum)');
%set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');


% L2/3, Model
subplot(1, 3, 2); hold on;
if (Parset_i ~= 3) & (Parset_i ~= 6), l = zeros(1, 4); lgdtxt = cell(1, 4);
else, l = zeros(1, 2); lgdtxt = cell(1, 2); end
Clr_L23 = [0.8 0 0.4; 0.5 0 0.8];
%
% Model, All pairs
l(1) = errorbar(d_BinCenter_E, NoiseCorr_mean_E, NoiseCorr_se_E,...
    'Color', Clr_L23(1, :), 'LineWidth', lwdth, 'Marker', '.',...
    'MarkerSize', mksize, 'CapSize', cpsize);
lgdtxt{1} = 'All pairs';
%
if (Parset_i ~= 3) & (Parset_i ~= 6)
    par = fitpar_noi_E; par_err = fitpar_err_noi_E;
    l(2) = plot(x2, DualExp2(x2, par(1), par(2), par(3), par(4), par(5)),...
        'LineStyle', '--', 'LineWidth', lwdth, 'Color', Clr_L23(1, :));
    lgdtxt{2} = ['(An = ', num2str(par(1), '%.2f'), ', \lambdan = ',...
        num2str(par(3), '%.2f'), ' \pm ', num2str(par_err(3), '%.2f'), ' \mum)'];
    %
    % Model, dtheta ~ 0
    l(3) = errorbar(d_BinCenter_E, NoiseCorr_mean_E_dtheta1, NoiseCorr_se_E_dtheta1,...
        'Color', Clr_L23(2, :), 'LineWidth', lwdth, 'Marker', '.',...
        'MarkerSize', mksize, 'CapSize', cpsize);
    lgdtxt{3} = ['\Delta\theta <= ', num2str(dtheta_bound1(2) * (180 / pi)), '^o'];
    %
    par = fitpar_noi_E_dtheta1; par_err = fitpar_err_noi_E_dtheta1;
    l(4) = plot(x2, DualExp2(x2, par(1), par(2), par(3), par(4), par(5)),...
        'LineStyle', '--', 'LineWidth', lwdth, 'Color', Clr_L23(2, :));
    lgdtxt{4} = ['(An = ', num2str(par(1), '%.2f'), ', \lambdan = ',...
        num2str(par(3), '%.2f'), ' \pm ', num2str(par_err(3), '%.2f'), ' \mum)'];
else
    % Model, dtheta ~ 0
    l(2) = errorbar(d_BinCenter_E, NoiseCorr_mean_E_dtheta1, NoiseCorr_se_E_dtheta1,...
        'Color', Clr_L23(2, :), 'LineWidth', lwdth, 'Marker', '.',...
        'MarkerSize', mksize, 'CapSize', cpsize);
    lgdtxt{2} = ['\Delta\theta <= ', num2str(dtheta_bound1(2) * (180 / pi)), '^o'];
end
%
plot([0, xylim(2)], [0 0], 'k--');
legend(l, lgdtxt, 'FontSize', txtsz);
title({'Layer 2/3 Model', ['(Time bin ', num2str(T_on), ' ms)']}, 'FontWeight', 'normal');
%
%axis(xylim); set(gca, 'XTick', xtickc, 'YTick', ytickc);
%axis square; grid on;
%xlabel('Horizontal Cortical Distance (\mum)');
%set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');


% L2/3, Data
subplot(1, 3, 3); hold on;
l = zeros(1, 4); lgdtxt = cell(1, 4);
Clr_L23 = [1 0 0; 0 0 1];
%
% Data, All pairs
par = FitPar(1, :); par_err = FitPar_err(1, :);
Corr_broad_est = DualExp2(d_BinCenter, 0, par(2), 114514, par(4), par(5));
Corr_narrow_est = Corr_avg{1} - Corr_broad_est;
l(1) = errorbar(d_BinCenter, Corr_narrow_est, Corr_se{1}, 'Color', Clr_L23(1, :),...
    'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
lgdtxt{1} = 'All pairs';
%
l(2) = plot(x2, DualExp2(x2, par(1), 0, par(3), 114514, 0),...
    'LineStyle', '--', 'LineWidth', lwdth, 'Color', Clr_L23(1, :));
lgdtxt{2} =  ['(An = ', num2str(par(1), '%.2f'), ', \lambdan = ',...
    num2str(par(3), '%.2f'), ' \pm ', num2str(par_err(3), '%.2f'), ' \mum)'];
%
% Data, dtheta ~ 0
par = FitPar(2, :); par_err = FitPar_err(2, :);
Corr_broad_est = DualExp2(d_BinCenter, 0, par(2), 114514, par(4), par(5));
Corr_narrow_est = Corr_avg{2} - Corr_broad_est;
l(3) = errorbar(d_BinCenter, Corr_narrow_est, Corr_se{2}, 'Color', Clr_L23(2, :),...
    'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
lgdtxt{3} = '\Delta\theta ~ same group';
%
l(4) = plot(x2, DualExp2(x2, par(1), 0, par(3), 114514, 0),...
    'LineStyle', '--', 'LineWidth', lwdth, 'Color', Clr_L23(2, :));
lgdtxt{4} =  ['(An = ', num2str(par(1), '%.2f'), ', \lambdan = ',...
    num2str(par(3), '%.2f'), ' \pm ', num2str(par_err(3), '%.2f'), ' \mum)'];
%
plot([0, xylim(2)], [0 0], 'k--');
legend(l, lgdtxt, 'FontSize', txtsz);
title({'Layer 2/3 Data', '(Narrow component only)'}, 'FontWeight', 'normal');
%
for k = 1: 3
subplot(1, 3, k);
axis(xylim); set(gca, 'XTick', xtickc, 'YTick', ytickc);
axis square; grid on;
xlabel('Horizontal Cortical Distance (\mum)');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
end
%
% yaxis of model
if (Parset_i == 1) | (Parset_i == 2) | (Parset_i == 3)
    subplot(1, 3, 2); ylim([-0.02 * (3 / 5) 0.3]); set(gca, 'YTick', 0: 0.025: 0.3);
elseif (Parset_i == 4) | (Parset_i == 5) | (Parset_i == 6)
    subplot(1, 3, 2); ylim([-0.02 * (3.5 / 5) 0.35]); set(gca, 'YTick', 0: 0.05: 0.35);
end
subplot(1, 3, 1); ylim([-0.02 * (10 / 5) 1]); set(gca, 'YTick', 0: 0.1: 1);


pause(2);
print(gcf, '-dpng', [dir_local, '/Figures/Results_C3_', num2str(Parset_i), '.png']);
close;



end

