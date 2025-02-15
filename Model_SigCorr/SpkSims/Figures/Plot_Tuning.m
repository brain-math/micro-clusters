Parset_i = [1, 2];
fitting_exp_type = 2;


%% Spatial model
%
load([pwd, '/Results/Results_Total_', num2str(Parset_i(1)), '.mat'], 'Pref_theta_F', 'd_tot',...
    'Pref_theta_E', 'Pref_theta_I', 'd_BinCenter_E', 'SigCorr_mean_E', 'SigCorr_se_E',...
    'd_BinCenter_I', 'SigCorr_mean_I', 'SigCorr_se_I', 'gOSI_E', 'gOSI_I');
if fitting_exp_type == 1
    load([pwd, '/Results/Results_Total_', num2str(Parset_i(1)), '.mat'],...
        'fitpar_exp1_E', 'CI_exp1_E');%, 'fitpar_exp1_I', 'CI_exp1_I');
    par_E = fitpar_exp1_E; par_E_err = (CI_exp1_E(:, 2) - CI_exp1_E(:, 1)) / 2;
    %par_I = fitpar_exp1_I; par_I_err = (CI_exp1_I(:, 2) - CI_exp1_I(:, 1)) / 2;
    Exp = @(x, A, lambda, b) A * exp(- x / lambda) + b; x0 = 7.5: 0.1: 50;
elseif fitting_exp_type == 2
    load([pwd, '/Results/Results_Total_', num2str(Parset_i(1)), '.mat'],...
        'fitpar_exp2_E', 'CI_exp2_E');%, 'fitpar_exp2_I', 'CI_exp2_I');
    par_E = fitpar_exp2_E; par_E_err = (CI_exp2_E(:, 2) - CI_exp2_E(:, 1)) / 2;
    %par_I = fitpar_exp2_I; par_I_err = (CI_exp2_I(:, 2) - CI_exp2_I(:, 1)) / 2;
    Exp = @(x, A, lambda, b) A * exp(- (x / lambda).^2) + b; x0 = 0: 0.1: 50;
end
%
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 0.72]);
%
%
Clb_hsv = hsv;
subplot(2, 5, 1); imagesc(Pref_theta_E); colormap(gca, Clb_hsv); axis square;
set(gca, 'XTick', [1 100], 'XTickLabel', {'0', [num2str(d_tot), ' \mum']}, 'YTick', [1 100]);
title('Pref. orientation, L2/3 E, spatial', 'fontweight', 'normal');
clb = colorbar; clb.Position = [0.09, 0.595, 0.01, 0.32];
set(clb, 'YTick', [0.002, 1/4: 1/4: 3/4, 1 - 0.002] * pi,...
    'YTickLabel', {'0^o', '45^o', '90^o', '135^o', '180^o'});
subplot(2, 5, 6); imagesc(Pref_theta_I); colormap(gca, Clb_hsv); axis square;
set(gca, 'XTick', [1 50], 'XTickLabel', {'0', [num2str(d_tot), ' \mum']}, 'YTick', [1 50]);
title('Pref. orientation, L2/3 I, spatial', 'fontweight', 'normal');
%
%
subplot(2, 5, 3); hold on;
l3 = zeros(1, 2);
l3(1) = errorbar(d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E, 'color', [1 0.5 0.5], 'linewidth', 1);
plot(x0, Exp(x0, par_E(1), par_E(2), par_E(3)), 'color', [1 0.5 0.5], 'lineStyle', '--', 'linewidth', 1);
%
subplot(2, 5, 8); hold on;
l8 = zeros(1, 2);
l8(1) = errorbar(d_BinCenter_I, SigCorr_mean_I, SigCorr_se_I, 'color', [0.5 0.5 1], 'linewidth', 1);
%plot(x0, Exp(x0, par_I(1), par_I(2), par_I(3)), 'color', [0.5 0.5 1], 'lineStyle', '--', 'linewidth', 1);
%
%
subplot(2, 5, 4); hold on;
lgdtxt4 = cell(1, 2);
gOSI_min = 0; gOSI_max = 1; d_gOSI = 0.02;
gOSI_ctr = gOSI_min: d_gOSI: gOSI_max;
gOSI_edge = (gOSI_min - d_gOSI/2): d_gOSI: (gOSI_max + d_gOSI/2);
gOSI_E_pdf = histcounts(gOSI_E, gOSI_edge, 'Normalization', 'pdf');
bar(gOSI_ctr, gOSI_E_pdf, 1.01, 'EdgeColor', [0.75 0 0], 'FaceColor', [1 0.5 0.5]);
lgdtxt4{1} = ['Spatial, mean ', num2str(mean(gOSI_E), '%.3f')];
%
subplot(2, 5, 9); hold on;
lgdtxt9 = cell(1, 2);
gOSI_I_pdf = histcounts(gOSI_I, gOSI_edge, 'Normalization', 'pdf');
bar(gOSI_ctr, gOSI_I_pdf, 1.01, 'EdgeColor', [0 0 0.75], 'FaceColor', [0.5 0.5 1]);
lgdtxt9{1} = ['Spatial, mean ', num2str(mean(gOSI_I), '%.3f')];
%
%
subplot(2, 5, 5); hold on;
ptheta_min = 0; ptheta_max = 180; d_ptheta = 7.5;
ptheta_ctr = (ptheta_min + d_ptheta/2): d_ptheta: (ptheta_max - d_ptheta/2);
ptheta_edge = ptheta_min: d_ptheta: ptheta_max;
ptheta_E_pdf = histcounts(Pref_theta_E(:) * 180 / pi, ptheta_edge, 'Normalization', 'pdf');
bar(ptheta_ctr, ptheta_E_pdf, 0.9, 'EdgeColor', [0.75 0 0], 'FaceColor', [1 0.5 0.5]);
%
subplot(2, 5, 10); hold on;
ptheta_I_pdf = histcounts(Pref_theta_I(:) * 180 / pi, ptheta_edge, 'Normalization', 'pdf');
bar(ptheta_ctr, ptheta_I_pdf, 0.9, 'EdgeColor', [0 0 0.75], 'FaceColor', [0.5 0.5 1]);







%% Ftr. Spc. model
%
load([pwd, '/Results/Results_Total_', num2str(Parset_i(2)), '.mat'],...
    'Pref_theta_E', 'Pref_theta_I', 'd_BinCenter_E', 'SigCorr_mean_E', 'SigCorr_se_E',...
    'd_BinCenter_I', 'SigCorr_mean_I', 'SigCorr_se_I', 'gOSI_E', 'gOSI_I');
if fitting_exp_type == 1
    load([pwd, '/Results/Results_Total_', num2str(Parset_i(2)), '.mat'],...
        'fitpar_exp1_E', 'CI_exp1_E');%, 'fitpar_exp1_I', 'CI_exp1_I');
    par_E = fitpar_exp1_E; par_E_err = (CI_exp1_E(:, 2) - CI_exp1_E(:, 1)) / 2;
    %par_I = fitpar_exp1_I; par_I_err = (CI_exp1_I(:, 2) - CI_exp1_I(:, 1)) / 2;
    Exp = @(x, A, lambda, b) A * exp(- x / lambda) + b; x0 = 7.5: 0.1: 50;
elseif fitting_exp_type == 2
    load([pwd, '/Results/Results_Total_', num2str(Parset_i(2)), '.mat'],...
        'fitpar_exp2_E', 'CI_exp2_E');%, 'fitpar_exp2_I', 'CI_exp2_I');
    par_E = fitpar_exp2_E; par_E_err = (CI_exp2_E(:, 2) - CI_exp2_E(:, 1)) / 2;
    %par_I = fitpar_exp2_I; par_I_err = (CI_exp2_I(:, 2) - CI_exp2_I(:, 1)) / 2;
    Exp = @(x, A, lambda, b) A * exp(- (x / lambda).^2) + b; x0 = 0: 0.1: 50;
end
%
%
subplot(2, 5, 2); imagesc(Pref_theta_E); colormap(gca, Clb_hsv); axis square;
set(gca, 'XTick', [1 100], 'XTickLabel', {'0', [num2str(d_tot), ' \mum']}, 'YTick', [1 100]);
title('Pref. orientation, L2/3 E, ftr. spc.', 'fontweight', 'normal');
subplot(2, 5, 7); imagesc(Pref_theta_I); colormap(gca, Clb_hsv); axis square;
set(gca, 'XTick', [1 50], 'XTickLabel', {'0', [num2str(d_tot), ' \mum']}, 'YTick', [1 50]);
title('Pref. orientation, L2/3 I, ftr. spc.', 'fontweight', 'normal');
%
%
subplot(2, 5, 3); hold on;
l3(2) = errorbar(d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E, 'color', [1 0 0], 'linewidth', 1);
plot(x0, Exp(x0, par_E(1), par_E(2), par_E(3)), 'color', [1 0 0], 'lineStyle', '--', 'linewidth', 1);
plot([0, 50], [0 0], 'k--');
axis square; grid on;
axis([0, 50, -0.1, 0.5]); set(gca, 'XTick', 0: 5: 50, 'YTick', -0.1: 0.1: 0.5);
legend(l3, {'Spatial', 'Ftr. spc.'});
title('Signal corr., L2/3 E', 'fontweight', 'normal');
%
subplot(2, 5, 8); hold on;
l8(2) = errorbar(d_BinCenter_I, SigCorr_mean_I, SigCorr_se_I, 'color', [0 0 1], 'linewidth', 1);
%plot(x0, Exp(x0, par_I(1), par_I(2), par_I(3)), 'color', [0 0 1], 'lineStyle', '--', 'linewidth', 1);
plot([0, 50], [0 0], 'k--');
axis square; grid on;
axis([0, 50, -0.05, 0.15]); set(gca, 'XTick', 0: 5: 50, 'YTick', -0.05: 0.025: 0.15);
legend(l8, {'Spatial', 'Ftr. spc.'});
xlabel('Horizontal distance (\mum)'); title('Signal corr., L2/3 I', 'fontweight', 'normal');
%
%
subplot(2, 5, 4); hold on;
gOSI_E_pdf = histcounts(gOSI_E, gOSI_edge, 'Normalization', 'pdf');
bar(gOSI_ctr, gOSI_E_pdf, 1.01, 'EdgeColor', [0.75 0 0], 'FaceColor', [1 0 0]);
lgdtxt4{2} = ['Ftr spc., mean ', num2str(mean(gOSI_E), '%.3f')]; legend(lgdtxt4);
axis square; grid on; axis([0 1 0 4]); set(gca, 'XTick', 0: 0.2: 1);
title('p.d.f. of gOSI, L2/3 E', 'fontweight', 'normal');
%
subplot(2, 5, 9); hold on;
gOSI_I_pdf = histcounts(gOSI_I, gOSI_edge, 'Normalization', 'pdf');
bar(gOSI_ctr, gOSI_I_pdf, 1.01, 'EdgeColor', [0 0 0.75], 'FaceColor', [0 0 1]);
lgdtxt9{2} = ['Ftr spc., mean ', num2str(mean(gOSI_I), '%.3f')]; legend(lgdtxt9);
axis square; grid on; axis([0 1 0 8]); set(gca, 'XTick', 0: 0.2: 1);
title('p.d.f. of gOSI, L2/3 I', 'fontweight', 'normal');
%
%
subplot(2, 5, 5); hold on;
ptheta_E_pdf = histcounts(Pref_theta_E(:) * 180 / pi, ptheta_edge, 'Normalization', 'pdf');
bar(ptheta_ctr, ptheta_E_pdf, 0.9, 'EdgeColor', [0.75 0 0], 'FaceColor', [1 0 0]);
legend('Spatial', 'Ftr. spc.');
axis square; grid on; axis([0 180 0 0.01]); set(gca, 'XTick', 0: 30: 180, 'YTick', 0: 0.002: 0.01);
xlabel('(deg.)'); title('p.d.f. of \theta\_pref., L2/3 E', 'fontweight', 'normal');
%
subplot(2, 5, 10); hold on;
ptheta_I_pdf = histcounts(Pref_theta_I(:) * 180 / pi, ptheta_edge, 'Normalization', 'pdf');
bar(ptheta_ctr, ptheta_I_pdf, 0.9, 'EdgeColor', [0 0 0.75], 'FaceColor', [0 0 1]);
legend('Spatial', 'Ftr. spc.');
axis square; grid on; axis([0 180 0 0.01]); set(gca, 'XTick', 0: 30: 180, 'YTick', 0: 0.002: 0.01);
xlabel('(deg.)'); title('p.d.f. of \theta\_pref., L2/3 I', 'fontweight', 'normal');


pause(2); print(gcf, '-dpng', [pwd, '/Figures/Tuning_',...
    num2str(Parset_i(1)), '_', num2str(Parset_i(2)), '.png']); %close;
print(gcf, '-depsc2', '-r1500', [pwd, '/Figures/Tuning_',...
    num2str(Parset_i(1)), '_', num2str(Parset_i(2)), '.eps']);











%load(savename_L4, 'd0_F', 'gOSI_F', 'Pref_theta_F', 'fitpar_exp2_F', 'CI_exp2_F', 'fitpar_exp1_F', 'CI_exp1_F');
%load(savename_result_after, 'gOSI_E', 'Pref_theta_E', 'gOSI_I', 'Pref_theta_I');
%if ifBroadOnly == 0
%load(savename_result_after, 'd_BinCenter_E', 'SigCorr_mean_E', 'SigCorr_se_E',...
%    'fitpar_exp2_E', 'CI_exp2_E', 'fitpar_exp1_E', 'CI_exp1_E',...
%    'd_BinCenter_I', 'SigCorr_mean_I', 'SigCorr_se_I',...
%    'fitpar_exp2_I', 'CI_exp2_I', 'fitpar_exp1_I', 'CI_exp1_I');
%end
%if ifBroadOnly == 0, xmax = 50; ymin = -0.1; dx = 5;
%elseif ifBroadOnly == 1, xmax = 450; ymin = -0.3; dx = 75; end
%%
%figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.075, 0, 0.55, 1]);
%subplot(3, 3, 1); histogram(gOSI_F, 'normalization', 'pdf'); xlim([0 1]); axis square;
%ylabel('L4'); title('gOSI', 'FontWeight', 'normal');
%subplot(3, 3, 2); imagesc(Pref_theta_F); colormap('hsv'); axis square;
%title(['\Deltax = ', num2str(d0_F), ' (\mum)'], 'FontWeight', 'normal');
%subplot(3, 3, 3); hold on;
%plot([0 xmax], [0 0], 'k--');
%errorbar(d_BinCenter_F, SigCorr_mean_F, SigCorr_se_F, 'b');
%if fitting_exp_type == 1
%    par_F = fitpar_exp1_F; lambda_err = (CI_exp1_F(2, 2) - CI_exp1_F(2, 1)) / 2;
%    plot(x1, Exp1(x1, par_F(1), par_F(2), par_F(3)), 'b--');
%elseif fitting_exp_type == 2
%    par_F = fitpar_exp2_F; lambda_err = (CI_exp2_F(2, 2) - CI_exp2_F(2, 1)) / 2;
%    plot(x2, Exp2(x2, par_F(1), par_F(2), par_F(3)), 'b--');
%end
%axis square; grid on; axis([0 50 -0.1 1.01]); set(gca, 'XTick', 0: 5: 50, 'YTick', -0.1: 0.1: 1);
%title(['\lambda (L4) = ', num2str(par_F(2), '%.2f'), ' \pm ',...
%    num2str(lambda_err, '%.2f'), ' (\mum)'], 'FontWeight', 'normal');
%%
%subplot(3, 3, 4); histogram(gOSI_E); xlim([0 1]); axis square; ylabel('L2/3 E');
%subplot(3, 3, 5); imagesc(Pref_theta_E); colormap('hsv'); axis square;
%subplot(3, 3, 6); hold on;
%plot([0 xmax], [0 0], 'k--');
%errorbar(d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E, 'r');
%if ifBroadOnly == 0
%if fitting_exp_type == 1
%    par_E = fitpar_exp1_E; lambda_err = (CI_exp1_E(2, 2) - CI_exp1_E(2, 1)) / 2;
%    plot(x1, Exp1(x1, par_E(1), par_E(2), par_E(3)), 'r--');
%elseif fitting_exp_type == 2
%    par_E = fitpar_exp2_E; lambda_err = (CI_exp2_E(2, 2) - CI_exp2_E(2, 1)) / 2;
%    plot(x2, Exp2(x2, par_E(1), par_E(2), par_E(3)), 'r--');
%end
%end
%axis square; grid on; axis([0 xmax ymin 1.01]); set(gca, 'XTick', 0: dx: xmax, 'YTick', ymin: 0.1: 1);
%if ifBroadOnly == 0, title(['\lambda (L2/3E) = ', num2str(par_E(2), '%.2f'), ' \pm ',...
%    num2str(lambda_err, '%.2f'), ' (\mum)'], 'FontWeight', 'normal'); end
%%
%subplot(3, 3, 7); histogram(gOSI_I); xlim([0 1]); axis square; ylabel('L2/3 I');
%subplot(3, 3, 8); imagesc(Pref_theta_I); colormap('hsv'); axis square;
%subplot(3, 3, 9); hold on;
%plot([0 xmax], [0 0], 'k--');
%errorbar(d_BinCenter_I, SigCorr_mean_I, SigCorr_se_I, 'r');
%if ifBroadOnly == 0
%if fitting_exp_type == 1
%    par_I = fitpar_exp1_I; lambda_err = (CI_exp1_I(2, 2) - CI_exp1_I(2, 1)) / 2;
%    plot(x1, Exp1(x1, par_I(1), par_I(2), par_I(3)), 'r--');
%elseif fitting_exp_type == 2
%    par_I = fitpar_exp2_I; lambda_err = (CI_exp2_I(2, 2) - CI_exp2_I(2, 1)) / 2;
%    plot(x2, Exp2(x2, par_I(1), par_I(2), par_I(3)), 'r--');
%end
%end
%axis square; grid on; axis([0 xmax ymin 1.01]); set(gca, 'XTick', 0: dx: xmax, 'YTick', ymin: 0.1: 1);
%if ifBroadOnly == 0, title(['\lambda (L2/3I) = ', num2str(par_I(2), '%.2f'), ' \pm ',...
%    num2str(lambda_err, '%.2f'), ' (\mum)'], 'FontWeight', 'normal'); end
%%
%pause(2); print(gcf, '-dpng', [savename_figure, '_supp.png']); close;

