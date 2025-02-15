%% Run for the first time
if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
addpath(genpath([dir_0, '/Study/CompNeuro/Projects/Functions_simul']));
dir_local = [dir_0, '/Study/CompNeuro/Projects/Micro-clustering'];
cd(dir_local); addpath(genpath([pwd, '/functions'])); cd([pwd, '/Model_EPerturb']);
if ~exist([pwd, '/Parameters'], 'dir'), mkdir([pwd, '/Parameters']); end
if ~exist([pwd, '/Figures'], 'dir'), mkdir([pwd, '/Figures']); end
if ~exist([pwd, '/Results'], 'dir'), mkdir([pwd, '/Results']); end
addpath(genpath(pwd));
rng('shuffle');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dist_Bin = 10;
%figure; hold on;

Parset_i = 4;
%for Parset_i = 1: 4
load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param');
N_trial_tot_valid = param.N_trial_tot(2);
Ne = param.Ne; Ne_sqrt = sqrt(Ne); d_e = param.d_tot / Ne_sqrt;
%
tmp = param.PertIdx_Tot(:, param.N_spk_pert: param.N_spk_pert: end);
tmp = tmp(param.N_simul(1) + 1: end, param.N_trial_burn + 1: end)';
PertIdx_ordered = tmp(:)'; clear tmp    % (1, N_trial_tot_valid)
%
d_circ = @(d, T) acos(cos(d * (2 * pi / T))) * (T / (2 * pi));
[jj, ii] = meshgrid([1: Ne_sqrt]); j_grid = jj(:); i_grid = ii(:); clear ii jj
r_to_PertIdx = NaN(Ne, N_trial_tot_valid);
for trial_k = 1: N_trial_tot_valid
    k_pert = PertIdx_ordered(trial_k);
    i_pert = mod(k_pert - 1, Ne_sqrt) + 1; j_pert = ceil(k_pert / Ne_sqrt);
    r_to_PertIdx(:, trial_k) = sqrt(d_circ(i_grid - i_pert, Ne_sqrt) .^ 2 +...
        d_circ(j_grid - j_pert, Ne_sqrt) .^ 2) * d_e;
    % imagesc(reshape(dist_k, [Ne_sqrt, Ne_sqrt])); colorbar;
end
%
if Parset_i ~= 4
    load([pwd, '/Results/Results_Total_', num2str(Parset_i), '.mat'], 'FR_E_on');    % (Ne, N_trial_tot_valid)
elseif Parset_i == 4
    load([pwd, '/Results/Results_Total_', num2str(Parset_i), '_A.mat'], 'FR_E_on_1');
    load([pwd, '/Results/Results_Total_', num2str(Parset_i), '_B.mat'], 'FR_E_on_2');
    FR_E_on = [FR_E_on_1, FR_E_on_2]; clear FR_E_on_1 FR_E_on_2
end
%
dist_max = ceil(max(r_to_PertIdx(:)) / dist_Bin) * dist_Bin;
dist_ctr = [dist_Bin: dist_Bin: dist_max]' + dist_Bin/2;
if Parset_i == 3
    dist_edge = [0, 12, 22, 30: dist_Bin: (dist_max + dist_Bin)];
    % [7.5], [10.6066], [15, 16.7705], [21.2132, 22.5, ...] => [7.5, 10.6066], [15, 16.7705, 21.2132], [22.5, ...]
    % Stupid way to solve why kappa0 still has a weak but sig. narrow peak.
elseif Parset_i ~= 3
    dist_edge = [0: dist_Bin: (dist_max + dist_Bin)];
end
%
[dFR_E_dist_mean, dFR_E_dist_se, ~] = histogram_mean_sem(FR_E_on(:), r_to_PertIdx(:), dist_edge);
dFR_E_dist_mean = dFR_E_dist_mean(2: end); dFR_E_dist_se = dFR_E_dist_se(2: end);
clear FR_E_on r_to_PertIdx
%
% load([pwd, '/Results/Results_Total_', num2str(Parset_i), '.mat'], 'FR_E_off');    % scalar
% dFR_E_dist_mean = dFR_E_dist_mean - FR_E_off;
dFR_E_dist_mean_0 = mean(dFR_E_dist_mean((dist_ctr >= 120) & (dist_ctr <= 400)));
dFR_E_dist_mean = dFR_E_dist_mean - dFR_E_dist_mean_0;
%
save([pwd, '/Results/Results_Total_', num2str(Parset_i), '_Ana.mat'],...
    'param', 'dist_ctr', 'dFR_E_dist_mean', 'dFR_E_dist_se');
%
figure;
errorbar(dist_ctr, dFR_E_dist_mean, dFR_E_dist_se); xlim([0 200]);
%
%end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Single Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
%
dist_Bin = 10;
ymin_left = -0.0045; ymax_left = 0.015; dy_left = 0.0015;
ymin_right = -0.03; ymax_right = 0.1; dy_right = 0.01;
%
Parset_i = 4;
%for Parset_i = [1 2 4]
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.45, 0.85]); hold on;
xmax = 200; l = zeros(1, 3);
%
yyaxis left;
load([pwd, '/Results/Results_Total_', num2str(Parset_i), '_Ana.mat'],...
    'param', 'dist_ctr', 'dFR_E_dist_mean', 'dFR_E_dist_se');
param1 = param; clear param
l(1) = errorbar(dist_ctr, dFR_E_dist_mean, dFR_E_dist_se,...
    'color', 'b', 'linewidth', 2);
%
plot([0 xmax], [0 0], 'k--');
xlim([0 xmax]); set(gca, 'xtick', 0: 10: xmax);
ylim([ymin_left ymax_left]); set(gca, 'ytick', ymin_left: dy_left: ymax_left);
xlabel('Distance to simulus position (\mum)'); ylabel('\DeltaFR (Hz)');
ax = gca; ax.YColor = 'b';
%
yyaxis right;
% Theory
% Based no dFR_E(d) = Wee(d) + (Wee * Wee)(d) + (Wei * Wie)(d) and lots of term drops
G = @(d, sigma) exp(-(d .^ 2) / (2 * sigma ^ 2)) / (2 * pi * sigma ^ 2);
dFR_E_theory = @(d, A, kappa_E, sigma_n_E, Wee, Wi, sigma_b_E, sigma_b_I)...
    A * (kappa_E * G(d, sigma_n_E) + G(d, sigma_b_E) +...
        Wee * G(d, sqrt(2) * sigma_b_E) - Wi * G(d, sqrt(sigma_b_E ^ 2 + sigma_b_I ^ 2)));
% Wi = Wei * Wie / Wee
d0 = 0: 0.1: xmax;
l(2) = plot(d0, dFR_E_theory(d0, 1400, param1.kappa(2), param1.sigma_micron(2, 2),...
    param1.W_est(1, 2), - param1.W_est(1, 3) * param1.W_est(2, 2) / param1.W_est(1, 2),...
    param1.sigma_micron(1, 2), param1.sigma_micron(1, 3)),...
    'color', [0.5 0 1], 'linewidth', 1, 'linestyle', '-', 'marker', 'none');
%
load([dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_EPerturb_OneCellATime',...
    '/OneCellAtATime_112N4M.mat'], 'x_d10', 'y_d10', 'y_se_d10');
x_data = x_d10; y_data = y_d10; y_se_data = y_se_d10;
l(3) = errorbar(x_data, y_data, y_se_data, 'color', 'r', 'linewidth', 2, 'linestyle', '-');
ylim([ymin_right ymax_right]); set(gca, 'ytick', ymin_right: dy_right: ymax_right);
ylabel('\DeltaF/F0');
xlabel('Distance to simulus position (\mum)'); ylabel('\DeltaF/F0');
ax = gca; ax.YColor = 'r';
axis square; grid on;
%
legend(l, {['Simulation (\kappaE = ', num2str(param1.kappa(2)),...
    ', ', num2str(param1.N_trial_tot(2)), ' trials)'],...
    'Theory', 'Data (112 neuron, 4 mice)'}, 'fontsize', 14);
set(gca, 'FontSize', 12);
%
s = param1.sigma_micron; k = param1.kappa; K = param1.K_in;
%
%sigma_text = ['\sigma = [144, 146, 110; 3.75, 20, 9.25] (\mum)'];
sigma_text = ['\sigma = [', num2str(s(1, 1)), ', ', num2str(s(1, 2)), ', ', num2str(s(1, 3)),...
    '; ', num2str(s(2, 1)), ', ', num2str(s(2, 2)), ', ', num2str(s(2, 3)), '] (\mum)'];
%kappa_text = ['\kappa = [0.075, 0.125, 0.125]'];
kappa_text = ['\kappa = [', num2str(k(1)), ', ', num2str(k(2)), ', ', num2str(k(3)), ']'];
%K_in_text = ['K\_in = [1000, 200, 200]'];
K_in_text = ['K\_in = [', num2str(K(1, 1)), ', ', num2str(K(1, 2)), ', ', num2str(K(1, 3)), ']'];
txt3 = [num2str(param1.alpha_slow * 100), '% slow PSP; ',...
    num2str(param1.T_off / 1000), 's off + ', num2str(param1.T_on / 1000),...
    's stimulation of ', num2str(param1.N_spk_pert / param1.T_on), ' Hz'];
%
title({sigma_text, [kappa_text, ', ', K_in_text], txt3}, 'FontWeight', 'normal');
%
savename = [pwd, '/Figures/Results_Par_', num2str(Parset_i)];%, '_d', num2str(dist_Bin)];
pause(2);
savefig([savename, '.fig']);
print(gcf, '-dpng', [savename, '.png']);
close;
%
%end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Group Plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
%
%Parset_i_kappa = 1; Parset_i_null = 3;
Parset_i_kappa = 4; Parset_i_null = 3;
%
dist_Bin = 10;
ymin_left = -0.0045; ymax_left = 0.015; dy_left = 0.0015;
ymin_right = -0.03; ymax_right = 0.1; dy_right = 0.01;
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.45, 0.85]); hold on;
xmax = 200; l = zeros(1, 4);
%
yyaxis left;
load([pwd, '/Results/Results_Total_', num2str(Parset_i_kappa), '_Ana.mat'],...
    'param', 'dist_ctr', 'dFR_E_dist_mean', 'dFR_E_dist_se');
param1 = param; clear param
l(1) = errorbar(dist_ctr, dFR_E_dist_mean, dFR_E_dist_se,...
    'color', 'b', 'linewidth', 2);
%
load([pwd, '/Results/Results_Total_', num2str(Parset_i_null), '_Ana.mat'],...
    'param', 'dist_ctr', 'dFR_E_dist_mean', 'dFR_E_dist_se');
param2 = param; clear param
l(2) = errorbar(dist_ctr, dFR_E_dist_mean, dFR_E_dist_se,...
    'color', 'b', 'linestyle', '-.', 'linewidth', 1);
%
plot([0 xmax], [0 0], 'k--');
xlim([0 xmax]); set(gca, 'xtick', 0: 10: xmax);
ylim([ymin_left ymax_left]); set(gca, 'ytick', ymin_left: dy_left: ymax_left);
xlabel('Distance to simulus position (\mum)'); ylabel('\DeltaFR (Hz)');
ax = gca; ax.YColor = 'b';
%
yyaxis right;
% Theory
% Based no dFR_E(d) = Wee(d) + (Wee * Wee)(d) + (Wei * Wie)(d) and lots of term drops
G = @(d, sigma) exp(-(d .^ 2) / (2 * sigma ^ 2)) / (2 * pi * sigma ^ 2);
dFR_E_theory = @(d, A, kappa_E, sigma_n_E, Wee, Wi, sigma_b_E, sigma_b_I)...
    A * (kappa_E * G(d, sigma_n_E) + G(d, sigma_b_E) +...
        Wee * G(d, sqrt(2) * sigma_b_E) - Wi * G(d, sqrt(sigma_b_E ^ 2 + sigma_b_I ^ 2)));
% Wi = Wei * Wie / Wee
d0 = 0: 0.1: xmax;
l(3) = plot(d0, dFR_E_theory(d0, 1400, param1.kappa(2), param1.sigma_micron(2, 2),...
    param1.W_est(1, 2), - param1.W_est(1, 3) * param1.W_est(2, 2) / param1.W_est(1, 2),...
    param1.sigma_micron(1, 2), param1.sigma_micron(1, 3)),...
    'color', [0.5 0 1], 'linewidth', 1, 'linestyle', '-', 'marker', 'none');
%
load([dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_EPerturb_OneCellATime',...
    '/OneCellAtATime_112N4M.mat'], 'x_d10', 'y_d10', 'y_se_d10');
x_data = x_d10; y_data = y_d10; y_se_data = y_se_d10;
l(4) = errorbar(x_data, y_data, y_se_data, 'color', 'r', 'linewidth', 2, 'linestyle', '-');
ylim([ymin_right ymax_right]); set(gca, 'ytick', ymin_right: dy_right: ymax_right);
ylabel('\DeltaF/F0');
xlabel('Distance to simulus position (\mum)'); ylabel('\DeltaF/F0');
ax = gca; ax.YColor = 'r';
axis square; grid on;
%
legend(l, {['Simulation (\kappaE = ', num2str(param1.kappa(2)),...
    ', ', num2str(param1.N_trial_tot(2)), ' trials)'],...
    ['Simulation (\kappa   = ', num2str(param2.kappa(2)),...
    ',        ', num2str(param2.N_trial_tot(2)), ' trials)'],...
    'Theory', 'Data (112 neuron, 4 mice)'}, 'fontsize', 14);
set(gca, 'FontSize', 12);
%
s = param1.sigma_micron; k = param1.kappa; K = param1.K_in;
%
%sigma_text = ['\sigma = [144, 146, 110; 3.75, 20, 9.25] (\mum)'];
sigma_text = ['\sigma = [', num2str(s(1, 1)), ', ', num2str(s(1, 2)), ', ', num2str(s(1, 3)),...
    '; ', num2str(s(2, 1)), ', ', num2str(s(2, 2)), ', ', num2str(s(2, 3)), '] (\mum)'];
%kappa_text = ['\kappa = [0.075, 0.125, 0.125]'];
kappa_text = ['\kappa = [', num2str(k(1)), ', ', num2str(k(2)), ', ', num2str(k(3)), ']'];
%K_in_text = ['K\_in = [1000, 200, 200]'];
K_in_text = ['K\_in = [', num2str(K(1, 1)), ', ', num2str(K(1, 2)), ', ', num2str(K(1, 3)), ']'];
txt3 = [num2str(param1.alpha_slow * 100), '% slow PSP; ',...
    num2str(param1.T_off / 1000), 's off + ', num2str(param1.T_on / 1000),...
    's stimulation of ', num2str(param1.N_spk_pert / param1.T_on), ' Hz'];
%
title({sigma_text, [kappa_text, ', ', K_in_text], txt3}, 'FontWeight', 'normal');
%
savename = [pwd, '/Figures/Results_Par_', num2str(Parset_i_kappa), '_', num2str(Parset_i_null)];
pause(2);
savefig([savename, '.fig']);
print(gcf, '-dpng', [savename, '.png']);
close;















%% Raster plot example
%load([pwd, '/Results/Results_Trial1.mat']);
%figure; scatter(spktrain(1, :), spktrain(2, :), 0.5, 'k', 'filled'); axis ij; %ylim([1 12500]);
%% 4 trials, 20000ms
%axis([0 20000 0 10000]);
%set(gca, 'xtick', [0 4000 5000 9000 10000 14000 15000 19000 20000]);


%Twin = 10; psth = zeros(1, 20000 / Twin);
%spk_time_bin_id = ceil(spktrain(1, :) / Twin);
%for k = 1: length(psth)
%    psth(k) = sum((spk_time_bin_id > k - 0.5) & (spk_time_bin_id < k + 0.5)) / (10000 * Twin / 1000);
%end
%figure; plot(Twin: Twin: 20000, psth)
%set(gca, 'xtick', [0 4000 5000 9000 10000 14000 15000 19000 20000]);







% % Theory 2: dFR = (eye(N) - W) \ (W * r_pert)    % r_pert = zeros(N, 1) except that excited neuron for 10Hz.
% % But that's poor theory based on ramdomly sampled W.
%load([pwd, '/Parameters/Parameters_Recurrent_1.mat'], 'param'); 
%Ne = param.Ne; Ne_sqrt = sqrt(Ne); Ni = param.Ni; N = param.N;
%d_tot = param.d_tot; r_pert_0 = param.N_spk_pert / param.T_on; clear param    % 10Hz
%Ge = 12; Gi = 15;
%G = diag([Ge * ones(1, Ne), Gi * ones(1, Ni)]);
%%
%Pert_i = 51; Pert_j = 51;
%[xx, yy] = meshgrid([1: Ne_sqrt] - Pert_j, [1: Ne_sqrt] - Pert_i);
%r_2d = sqrt(xx .^ 2 + yy .^ 2) * (d_tot / Ne_sqrt); clear xx yy
%%
%r_pert = zeros(N, 1); r_pert((Pert_j - 1) * Ne_sqrt + Pert_i) = r_pert_0;


%Parset_i = 1;
%load([pwd, '/Parameters/JN_', num2str(Parset_i), '.mat'], 'JN_ee', 'JN_ei', 'JN_ie', 'JN_ii');
%W = G * [JN_ee, JN_ei; JN_ie, JN_ii]; clear JN_ee JN_ei JN_ie JN_ii
%dFR_theory = (eye(N) - W) \ (W * r_pert);
%dFR_E_theory = dFR_theory(1: Ne); clear W dFR_theory

%dist_Bin = 10; dist_max = ceil(max(r_2d(:)) / dist_Bin) * dist_Bin;
%dist_ctr = [dist_Bin: dist_Bin: dist_max]' + dist_Bin/2;
%dist_edge = [0: dist_Bin: (dist_max + dist_Bin)];
%%
%[dFR_E_dist_mean_theory, dFR_E_dist_se_theory, ~] = histogram_mean_sem(dFR_E_theory(:), r_2d(:), dist_edge);
%dFR_E_dist_mean_theory = dFR_E_dist_mean_theory(2: end);
%dFR_E_dist_se_theory = dFR_E_dist_se_theory(2: end);
%%
%dFR_E_dist_mean_theory_0 = mean(dFR_E_dist_mean_theory((dist_ctr >= 200) & (dist_ctr <= 400)));
%dFR_E_dist_mean_theory = dFR_E_dist_mean_theory - dFR_E_dist_mean_theory_0;

%figure; errorbar(dist_ctr, dFR_E_dist_mean_theory, dFR_E_dist_se_theory); xlim([0 200]);

