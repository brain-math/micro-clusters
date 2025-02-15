function plotL4(d_BinCenter_F, SigCorr_mean_F, SigCorr_se_F,...
    fitpar_exp2_F, CI_exp2_F, fitpar_exp1_F, CI_exp1_F, Exp_type)

if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
load([dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging',...
    '/Analysis_5_JointStat_Overall/SigCorr_plot_d5.mat'],...
    'd_BinStep_overlap', 'd_BinCenter_data2', 'SigCorr_mean_data2', 'SigCorr_se_data2',...
    'fitpar_exp1_data2', 'fitpar_exp1_err_data2', 'x1', 'Exp1',...
    'fitpar_exp2_data2', 'fitpar_exp2_err_data2', 'x2', 'Exp2', 'Clb_2');

data_type_title = {'L2/3', 'L4'};
data_type_list = [2];
l = zeros(1, 2 * length(data_type_list) + 2);
lgdtxt = cell(1, 2 * length(data_type_list) + 2);
%
lwdth = 1.5; mksize = 12.5; cpsize = 7.5; txtsz = 15;
xylim = [0 50 -0.1 1.01]; xtickc = [0: 5: 50]; ytickc = -0.1: 0.1: 1;
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.55, 1]); hold on;
suptitle({'Signal Correlation ~ Horizontal Cortical distance',...
        ['(Bin Size = ', num2str(d_BinStep_overlap), ' \mum)']}, 8, 0.95);
for i = 1: length(data_type_list)
    k = data_type_list(i);
    l(2 * i - 1) = errorbar(d_BinCenter_data2, SigCorr_mean_data2{k}, SigCorr_se_data2{k},...
        'Color', Clb_2(k, :), 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
    lgdtxt{2 * i - 1} = data_type_title{k};
    if Exp_type == 1, par = fitpar_exp1_data2; par_err = fitpar_exp1_err_data2; 
        l(2 * i) = plot(x1, Exp1(x1, par(k, 1), par(k, 2), par(k, 3)), 'LineStyle', '-.', 'Color', Clb_2(k, :));
    elseif Exp_type == 2, par = fitpar_exp2_data2; par_err = fitpar_exp2_err_data2; 
        l(2 * i) = plot(x2, Exp2(x2, par(k, 1), par(k, 2), par(k, 3)), 'LineStyle', '-.', 'Color', Clb_2(k, :));
    end
    lgdtxt{2 * i} = ['(\lambda = ', num2str(par(k, 2), '%.2f'), ' \pm ', num2str(par_err(k, 2), '%.2f'),...
        ' \mum , A = ', num2str(par(k, 1), '%.3f'), ', b = ', num2str(par(k, 3), '%.3f'), ')'];
end
plot([0, xylim(2)], [0 0], 'k--');
axis(xylim); set(gca, 'XTick', xtickc, 'YTick', ytickc);
axis square; grid on; 
xlabel('Horizontal Cortical Distance (\mum)');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
Clr_L4_model = [0, 0.5, 1];
l(end - 1) = errorbar(d_BinCenter_F, SigCorr_mean_F, SigCorr_se_F, 'Color', Clr_L4_model, 'LineWidth', lwdth);
lgdtxt{end - 1} = 'L4 model';
if Exp_type == 1
    l(end) = plot(x1, Exp1(x1, fitpar_exp1_F(1), fitpar_exp1_F(2), fitpar_exp1_F(3)),...
        'Color',  Clr_L4_model, 'LineStyle', '--', 'LineWidth', lwdth);
    lgdtxt{end} = ['(\lambda = ', num2str(fitpar_exp1_F(2), '%.2f'), ' \pm ',...
        num2str((CI_exp1_F(2, 2) - CI_exp1_F(2, 1)) / 2, '%.2f'),...
        ' \mum , A = ', num2str(fitpar_exp1_F(1), '%.3f'), ', b = ', num2str(fitpar_exp1_F(3), '%.3f'), ')'];
elseif Exp_type == 2
    l(end) = plot(x2, Exp2(x2, fitpar_exp2_F(1), fitpar_exp2_F(2), fitpar_exp2_F(3)),...
        'Color', Clr_L4_model, 'LineStyle', '--', 'LineWidth', lwdth);
    lgdtxt{end} = ['(\lambda = ', num2str(fitpar_exp2_F(2), '%.2f'), ' \pm ',...
        num2str((CI_exp2_F(2, 2) - CI_exp2_F(2, 1)) / 2, '%.2f'),...
        ' \mum , A = ', num2str(fitpar_exp2_F(1), '%.3f'), ', b = ', num2str(fitpar_exp2_F(3), '%.3f'), ')'];
end
legend(l, lgdtxt, 'FontSize', txtsz);

