function plot_spatial(savename_L4, savename_result, ifBroadOnly, fitting_exp_type, savename_figure)

load(savename_L4, 'd_BinCenter_F', 'SigCorr_mean_F', 'SigCorr_se_F');
load(savename_result, 'd_BinCenter_E', 'SigCorr_mean_E', 'SigCorr_se_E');

if ifBroadOnly == 0    % Then compare with data.
if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
load([dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging/',...
    'Analysis_5_JointStat_Overall/SigCorr_plot_d5.mat'],...
    'd_BinCenter_data2', 'SigCorr_mean_data2', 'SigCorr_se_data2',...
    'fitpar_exp1_data2', 'fitpar_exp1_err_data2', 'Exp1',...
    'fitpar_exp2_data2', 'fitpar_exp2_err_data2', 'Exp2', 'Clb_2');
load(savename_result, 'fitpar_exp2_E', 'CI_exp2_E', 'fitpar_exp1_E', 'CI_exp1_E');
x1 = 7.5: 0.1: 50; x2 = 0: 0.1: 50;
%
% plot order: L4 data, L4 model, L2/3 data, L2/3 model
% Data
data_type_title = {'L2/3 data', 'L4 data'};
legend_order_data = [1 3]; legend_order_model = [2 4];
l = zeros(1, length(data_type_title) + 2); lgdtxt = cell(1, length(data_type_title) + 2);
lwdth = 1.5; mksize = 12.5; cpsize = 7.5; txtsz = 15;
xylim = [0 50 -0.1 1.01]; xtickc = [0: 5: 50]; ytickc = -0.1: 0.1: 1;
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.55, 1]); hold on;
for k = 1: length(data_type_title)
    m = legend_order_data(k);
    l(m) = errorbar(d_BinCenter_data2, SigCorr_mean_data2{k}, SigCorr_se_data2{k},...
        'Color', Clb_2(k, :), 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
    if fitting_exp_type == 1
        par = fitpar_exp1_data2; par_err = fitpar_exp1_err_data2;
        plot(x1, Exp1(x1, par(1), par(2), par(3)), 'LineStyle', '--', 'LineWidth', lwdth, 'Color', Clb_2(k, :));
    elseif fitting_exp_type == 2
        par = fitpar_exp2_data2; par_err = fitpar_exp2_err_data2;
        plot(x2, Exp2(x2, par(1), par(2), par(3)), 'LineStyle', '--', 'LineWidth', lwdth, 'Color', Clb_2(k, :));
    end
    lgdtxt{m} = [data_type_title{k}, ' (\lambda = ', num2str(par(k, 2), '%.2f'),...
        ' \pm ', num2str(par_err(k, 2), '%.2f'), ' \mum)'];
end
plot([0, xylim(2)], [0 0], 'k--');
axis(xylim); set(gca, 'XTick', xtickc, 'YTick', ytickc);
axis square; grid on;
xlabel('Horizontal Cortical Distance (\mum)');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
% Model
Clb_L4_model = [0, 0.5, 1]; Clb_L23_model = [0.5, 0, 0.5];
l(legend_order_model(1)) = errorbar(d_BinCenter_F, SigCorr_mean_F, SigCorr_se_F,...
    'Color', Clb_L4_model, 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
lgdtxt{legend_order_model(1)} = 'L4 model input';
%
l(legend_order_model(2)) = errorbar(d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E,...
    'Color', Clb_L23_model, 'LineWidth', lwdth + 1, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
if fitting_exp_type == 1, par = fitpar_exp1_E; lambda_err = (CI_exp1_E(2, 2) - CI_exp1_E(2, 1)) / 2;
    plot(x1, Exp1(x1, par(1), par(2), par(3)), 'LineStyle', '--', 'LineWidth', lwdth + 1, 'Color', Clb_L23_model);
elseif fitting_exp_type == 2, par = fitpar_exp2_E; lambda_err = (CI_exp2_E(2, 2) - CI_exp2_E(2, 1)) / 2;
    plot(x2, Exp2(x2, par(1), par(2), par(3)), 'LineStyle', '--', 'LineWidth', lwdth + 1, 'Color', Clb_L23_model);
end
lgdtxt{legend_order_model(2)} = ['L2/3 model output (\lambda = ',...
    num2str(par(2), '%.2f'), ' \pm ', num2str(lambda_err, '%.2f'), ' \mum)'];
legend(l, lgdtxt, 'FontSize', txtsz);
%
load(savename_result, 'sigma_micron', 'kappa', 'W', 'K_in');
sigma_b_text = ['\sigma\_b = [', num2str(sigma_micron(1, 1)),...
    ', ', num2str(sigma_micron(1, 2)), ', ', num2str(sigma_micron(1, 3)), '];'];
sigma_n_text = ['\sigma\_n = [', num2str(sigma_micron(2, 1)),...
    ', ', num2str(sigma_micron(2, 2)), ', ', num2str(sigma_micron(2, 3)), '] (\mum);'];
kappa_text = ['\kappa = [', num2str(kappa(1)),...
    ', ', num2str(kappa(2)), ', ', num2str(kappa(3)), '];'];
W_text = ['W = [', num2str(W(1, 1)), ', ', num2str(W(1, 2)), ', ', num2str(W(1, 3)), '; ',...
    num2str(W(2, 1)), ', ', num2str(W(2, 2)), ', ', num2str(W(2, 3)), '];'];
K_text = ['K\_in = [', num2str(K_in(1, 1)), ', ', num2str(K_in(1, 2)), ', ', num2str(K_in(1, 3)), '];'];
suptitle({[sigma_b_text, ' ', sigma_n_text, ' ', kappa_text], [W_text, ' ', K_text]}, 8, 0.95);
% suptitle({['\sigma\_b = [', num2str(sigma_micron(1, 1)), ', ', num2str(sigma_micron(1, 2)),...
%     ', ', num2str(sigma_micron(1, 3)), ']; \sigma\_n = [', num2str(sigma_micron(2, 1)),...
%     ', ', num2str(sigma_micron(2, 2)), ', ', num2str(sigma_micron(2, 3)), '] (\mum); \kappa = [',...
%     num2str(kappa(1)), ', ', num2str(kappa(2)), ', ', num2str(kappa(3)), ']'],...
%     ['W = [', num2str(W(1, 1)), ', ', num2str(W(1, 2)), ', ', num2str(W(1, 3)),...
%     '; ', num2str(W(2, 1)), ', ', num2str(W(2, 2)), ', ', num2str(W(2, 3)), ']; K\_in = [',...
%     num2str(K_in(1, 1)), ', ', num2str(K_in(1, 2)), ', ', num2str(K_in(1, 3)), ']']}, 8, 0.95);
pause(2); print(gcf, '-dpng', [savename_figure, '.png']); %close;
end    % if ifBroadOnly == 0



% load(savename_L4, 'd0_F', 'gOSI_F', 'Pref_theta_F', 'fitpar_exp2_F', 'CI_exp2_F', 'fitpar_exp1_F', 'CI_exp1_F');
% load(savename_result, 'gOSI_E', 'Pref_theta_E', 'gOSI_I', 'Pref_theta_I');
% if ifBroadOnly == 0, load(savename_result, 'd_BinCenter_I', 'SigCorr_mean_I', 'SigCorr_se_I',...
%     'fitpar_exp2_I', 'CI_exp2_I', 'fitpar_exp1_I', 'CI_exp1_I'); end
% if ifBroadOnly == 0, xmax = 50; ymin = -0.1; dx = 5;
% elseif ifBroadOnly == 1, xmax = 450; ymin = -0.3; dx = 75; end
% %
% figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.075, 0, 0.55, 1]);
% subplot(3, 3, 1); histogram(gOSI_F, 'normalization', 'pdf'); xlim([0 1]); axis square;
% ylabel('L4'); title('gOSI', 'FontWeight', 'normal');
% subplot(3, 3, 2); imagesc(Pref_theta_F); colormap('hsv'); axis square;
% title(['\Deltax = ', num2str(d0_F), ' (\mum)'], 'FontWeight', 'normal');
% subplot(3, 3, 3); hold on;
% plot([0 xmax], [0 0], 'k--');
% errorbar(d_BinCenter_F, SigCorr_mean_F, SigCorr_se_F, 'b');
% if fitting_exp_type == 1
%     par_F = fitpar_exp1_F; lambda_err = (CI_exp1_F(2, 2) - CI_exp1_F(2, 1)) / 2;
%     plot(x1, Exp1(x1, par_F(1), par_F(2), par_F(3)), 'b--');
% elseif fitting_exp_type == 2
%     par_F = fitpar_exp2_F; lambda_err = (CI_exp2_F(2, 2) - CI_exp2_F(2, 1)) / 2;
%     plot(x2, Exp2(x1, par_F(1), par_F(2), par_F(3)), 'b--');
% end
% axis square; grid on; axis([0 50 -0.1 1.01]); set(gca, 'XTick', 0: 5: 50, 'YTick', -0.1: 0.1: 1);
% title(['\lambda (L4) = ', num2str(par_F(2), '%.2f'), ' \pm ',...
%     num2str(lambda_err, '%.2f'), ' (\mum)'], 'FontWeight', 'normal');
% %
% subplot(3, 3, 4); histogram(gOSI_E); xlim([0 1]); axis square; ylabel('L2/3 E');
% subplot(3, 3, 5); imagesc(Pref_theta_E); colormap('hsv'); axis square;
% subplot(3, 3, 6); hold on;
% plot([0 xmax], [0 0], 'k--');
% errorbar(d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E, 'r');
% if ifBroadOnly == 0
% if fitting_exp_type == 1
%     par_E = fitpar_exp1_E; lambda_err = (CI_exp1_E(2, 2) - CI_exp1_E(2, 1)) / 2;
%     plot(x1, Exp1(x1, par_E(1), par_E(2), par_E(3)), 'r--');
% elseif fitting_exp_type == 2
%     par_E = fitpar_exp2_E; lambda_err = (CI_exp2_E(2, 2) - CI_exp2_E(2, 1)) / 2;
%     plot(x2, Exp2(x2, par_E(1), par_E(2), par_E(3)), 'r--');
% end
% end
% axis square; grid on; axis([0 xmax ymin 1.01]); set(gca, 'XTick', 0: dx: xmax, 'YTick', ymin: 0.1: 1);
% if ifBroadOnly == 0, title(['\lambda (L2/3E) = ', num2str(par_E(2), '%.2f'), ' \pm ',...
%     num2str(lambda_err, '%.2f'), ' (\mum)'], 'FontWeight', 'normal'); end
% %
% subplot(3, 3, 7); histogram(gOSI_I); xlim([0 1]); axis square; ylabel('L2/3 I');
% subplot(3, 3, 8); imagesc(Pref_theta_I); colormap('hsv'); axis square;
% subplot(3, 3, 9); hold on;
% plot([0 xmax], [0 0], 'k--');
% errorbar(d_BinCenter_I, SigCorr_mean_I, SigCorr_se_I, 'r');
% if ifBroadOnly == 0
% if fitting_exp_type == 1
%     par_I = fitpar_exp1_I; lambda_err = (CI_exp1_I(2, 2) - CI_exp1_I(2, 1)) / 2;
%     plot(x1, Exp1(x1, par_I(1), par_I(2), par_I(3)), 'r--');
% elseif fitting_exp_type == 2
%     par_I = fitpar_exp2_I; lambda_err = (CI_exp2_I(2, 2) - CI_exp2_I(2, 1)) / 2;
%     plot(x2, Exp2(x2, par_I(1), par_I(2), par_I(3)), 'r--');
% end
% end
% axis square; grid on; axis([0 xmax ymin 1.01]); set(gca, 'XTick', 0: dx: xmax, 'YTick', ymin: 0.1: 1);
% if ifBroadOnly == 0, title(['\lambda (L2/3I) = ', num2str(par_I(2), '%.2f'), ' \pm ',...
%     num2str(lambda_err, '%.2f'), ' (\mum)'], 'FontWeight', 'normal'); end
% %
% pause(2); print(gcf, '-dpng', [savename_figure, '_supp.png']); close;

