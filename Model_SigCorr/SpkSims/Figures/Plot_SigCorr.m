%open([pwd, '/Figures/data_2.5.fig']);
%% click
%load([pwd, '/Results/Results_Total_3.mat']);
%errorbar(d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E);

%[mean(mean(L23_FR(1: Ne, :))), mean(mean(L23_FR(Ne + 1: end, :)))]



Parset_i = [1, 2]; d_data = 5;
Nf_sqrt = 100; fitting_exp_type = 2;



if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
savename_result = @(k) [dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Model_SigCorr',...
    '/SpkSims/Results/Results_Total_', num2str(Parset_i(k)), '.mat'];

load([dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging/',...
    'Analysis_5_JointStat_Overall/SigCorr_plot_d', num2str(d_data), '.mat'],...
    'd_BinCenter_data2', 'SigCorr_mean_data2', 'SigCorr_se_data2',...
    'fitpar_exp1_data2', 'fitpar_exp1_err_data2', 'Exp1',...
    'fitpar_exp2_data2', 'fitpar_exp2_err_data2', 'Exp2', 'Clb_2');
x1 = 7.5: 0.1: 50; x2 = 0: 0.1: 50;
%
% plot order: L4 data, L2/3 data, L4 model, L2/3 model (spatial), L2/3 model (feature spc)
N_legend = 5; legend_order_data = [1 2]; legend_order_model = [3 4 5];
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
load(savename_result(1), 'd_BinCenter_F', 'SigCorr_mean_F', 'SigCorr_se_F');
Clb_model = [0, 0.5, 1; 0.4, 0.4, 0.4; 0.5, 0, 0.5];
l(legend_order_model(1)) = errorbar(d_BinCenter_F, SigCorr_mean_F, SigCorr_se_F,...
    'Color', Clb_model(1, :), 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
lgdtxt{legend_order_model(1)} = 'L4 model input';
%
ldx_name_head = {'(spatial)', '(ftr. spc.)'};
for k = 1: 2
load(savename_result(k), 'Ne', 'L23_FR', 'd_BinCenter_E', 'SigCorr_mean_E', 'SigCorr_se_E',...
    'fitpar_exp2_E', 'CI_exp2_E', 'fitpar_exp1_E', 'CI_exp1_E');
r_pop_mean(k, :) = [mean(mean(L23_FR(1: Ne, :))), mean(mean(L23_FR(Ne + 1: end, :)))]; clear L23_FR
%
l(legend_order_model(k + 1)) = errorbar(d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E,...
    'Color', Clb_model(k + 1, :), 'LineWidth', lwdth + k - 1,...
    'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
if fitting_exp_type == 1
    par = fitpar_exp1_E; lambda_err = (CI_exp1_E(2, 2) - CI_exp1_E(2, 1)) / 2;
    plot(x1, Exp1(x1, par(1), par(2), par(3)), 'LineStyle', '--',...
        'LineWidth', lwdth + k - 1, 'Color', Clb_model(k + 1, :));
elseif fitting_exp_type == 2
    par = fitpar_exp2_E; lambda_err = (CI_exp2_E(2, 2) - CI_exp2_E(2, 1)) / 2;
    plot(x2, Exp2(x2, par(1), par(2), par(3)), 'LineStyle', '--',...
        'LineWidth', lwdth + k - 1, 'Color', Clb_model(k + 1, :));
end
lgdtxt{legend_order_model(k + 1)} = ['L2/3 model ', ldx_name_head{k}, ' (A = ',...
    num2str(par(1), '%.2f'), ', \lambda = ', num2str(par(2), '%.2f'),...
    ' \pm ', num2str(lambda_err, '%.2f'), ' \mum)'];
end
%
legend(l, lgdtxt, 'FontSize', txtsz);
%
load(savename_result(2), 'sigma_micron', 'kappa', 'K_in', 'ifDthetaDep', 'sigma_dtheta', 'a_dtheta', 'alpha_slow');
sigma_b_text = ['\sigmaB = [', num2str(sigma_micron(1, 1)),...
    ', ', num2str(sigma_micron(1, 2)), ', ', num2str(sigma_micron(1, 3)), '];'];
sigma_n_text = ['\sigmaN = [', num2str(sigma_micron(2, 1)),...
    ', ', num2str(sigma_micron(2, 2)), ', ', num2str(sigma_micron(2, 3)), '] (\mum);'];
kappa_text = ['\kappa = [', num2str(kappa(1)),...
    ', ', num2str(kappa(2)), ', ', num2str(kappa(3)), '];'];
K_text = ['K\_in = [', num2str(K_in(1, 1)), ', ', num2str(K_in(1, 2)), ', ', num2str(K_in(1, 3)), '];'];
alpha_text = ['\alpha\_slow = ', num2str(alpha_slow)];
r_text = ['r\_avg. = [', num2str(r_pop_mean(1, 1), '%.2f'),...
    ', ', num2str(r_pop_mean(1, 2), '%.2f'), '] -> [', num2str(r_pop_mean(2, 1), '%.2f'),...
    ', ', num2str(r_pop_mean(2, 2), '%.2f'), '] (Hz);'];
rewiring_text = ['Rewiring rule: [', num2str(ifDthetaDep(1, 1)), ', ', num2str(ifDthetaDep(1, 2)),...
    ', ', num2str(ifDthetaDep(1, 3)), '; ', num2str(ifDthetaDep(2, 1)), ', ', num2str(ifDthetaDep(2, 2)),...
    ', ', num2str(ifDthetaDep(2, 3)), '], \sigma = ', num2str(sigma_dtheta * (180 / pi)),...
    ' (deg), a = ', num2str(a_dtheta)];
suptitle({[sigma_b_text, ' ', sigma_n_text],...
    [kappa_text, ' ', K_text, ' ', alpha_text], rewiring_text, r_text}, 3, 0.95);


pause(2); print(gcf, '-dpng', [pwd, '/Figures/SigCorr_',...
    num2str(Parset_i(1)), '_', num2str(Parset_i(2)), '_d', num2str(d_data), '.png']);
savefig([pwd, '/Figures/SigCorr_',...
    num2str(Parset_i(1)), '_', num2str(Parset_i(2)), '_d', num2str(d_data), '.fig']);
%close;


