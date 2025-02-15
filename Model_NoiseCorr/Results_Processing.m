function Results_Processing

Parset_i = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%
load([pwd, '/Parameters/Parameters_FFWD.mat'], 'Nf', 'Nf_sqrt', 'd_tot', 'N_stim', 'theta_stim');
load([pwd, '/Parameters/Parameters_Recurrent_1.mat'], 'param');

% Each simulation has 50 valid trials.
% For each parameter,
  % 1 - 5: 45 degree stimulus in L4 for trial var, totally 250 trials.
  % 6 - 10: 90 degree.
  % 11 - 15: 135 degree.
  % 16 - 20: 180 degree. So totally 1000 trials for evoked.
%
N_Trial_per_simul = param.N_Trial_per_simul;    % 50
N_Simul = param.N_Simul_per_par;    % 20
N_trial_tot = N_Simul * N_Trial_per_simul;    % 1000
N_trial = N_trial_tot / N_stim;    % 250
%
Ne = param.Ne; Ne_sqrt = sqrt(Ne);
Ni = param.Ni; Ni_sqrt = sqrt(Ni); Nrec = param.N; clear param
%
d_BinStep_micron = 0;    % real ones


L4_Trials_evoked = zeros(Nf, N_trial, N_stim);
L23_Trials_evoked = zeros(Nrec, N_trial, N_stim);
%
N_repeat = N_Simul / N_stim;    % 5
for stim_k = 1: N_stim
for repeat_k = 1: N_repeat
    simul_k = (stim_k - 1) * N_repeat + repeat_k;
    load([pwd, '/Results/Results_Par', num2str(Parset_i), '_Simul',...
        num2str(simul_k), '.mat'], 'L4_FR', 'L23_FR');
    idx = (repeat_k - 1) * N_Trial_per_simul + [1: N_Trial_per_simul];
    L4_Trials_evoked(:, idx, stim_k) = L4_FR;
    L23_Trials_evoked(:, idx, stim_k) = L23_FR;
end
end
clear N_repeat stim_k repeat_k simul_k L4_FR L23_FR idx


% Signal correlation (as control test)
L4_FR_tuning = squeeze(mean(L4_Trials_evoked, 2));
L23E_FR_tuning = squeeze(mean(L23_Trials_evoked(1: Ne, :, :), 2));
L23I_FR_tuning = squeeze(mean(L23_Trials_evoked(Ne + 1: end, :, :), 2));
%
[gOSI_F, Pref_theta_F, d_BinCenter_F, SigCorr_mean_F, SigCorr_se_F, fitpar_sig_F, CI_sig_F, ~, ~] =...
    FR_Analysis(L4_FR_tuning, theta_stim, d_BinStep_micron, d_tot, Nf);
fitpar_err_sig_F = (CI_sig_F(:, 2) - CI_sig_F(:, 1)) / 2;
[gOSI_E, Pref_theta_E, d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E, fitpar_sig_E, CI_sig_E, ~, ~] =...
    FR_Analysis(L23E_FR_tuning, theta_stim, d_BinStep_micron, d_tot, Ne);
fitpar_err_sig_E = (CI_sig_E(:, 2) - CI_sig_E(:, 1)) / 2;
[gOSI_I, Pref_theta_I, d_BinCenter_I, SigCorr_mean_I, SigCorr_se_I, fitpar_sig_I, CI_sig_I, ~, ~] =...
    FR_Analysis(L23I_FR_tuning, theta_stim, d_BinStep_micron, d_tot, Ni);
fitpar_err_sig_I = (CI_sig_I(:, 2) - CI_sig_I(:, 1)) / 2;


% Noise correlation
[~, NoiseCorr_mean_F, NoiseCorr_se_F, fitpar_noi_F, CI_noi_F, ~]...
    = NoiseCorr_Analysis(L4_Trials_evoked, d_BinStep_micron, d_tot, 0);
fitpar_err_noi_F = (CI_noi_F(:, 2) - CI_noi_F(:, 1)) / 2;
%
[~, NoiseCorr_mean_E, NoiseCorr_se_E, fitpar_noi_E, CI_noi_E, ~]...
    = NoiseCorr_Analysis(L23_Trials_evoked(1: Ne, :, :), d_BinStep_micron, d_tot, 0);
fitpar_err_noi_E = (CI_noi_E(:, 2) - CI_noi_E(:, 1)) / 2;


% Pick up Δθ ~ 0 pairs for L2/3E
dtheta_bound1 = [0, 15] * (pi / 180);
[~, NoiseCorr_mean_E_dtheta1, NoiseCorr_se_E_dtheta1, fitpar_noi_E_dtheta1, CI_noi_E_dtheta1, ~] =...
    NoiseCorr_Analysis_dtheta(L23_Trials_evoked(1: Ne, :, :),...
    Pref_theta_E, dtheta_bound1, d_BinStep_micron, d_tot, 0);
fitpar_err_noi_E_dtheta1 = (CI_noi_E_dtheta1(:, 2) - CI_noi_E_dtheta1(:, 1)) / 2;
%
dtheta_bound2 = [0, 45] * (pi / 180);
[~, NoiseCorr_mean_E_dtheta2, NoiseCorr_se_E_dtheta2, fitpar_noi_E_dtheta2, CI_noi_E_dtheta2, ~] =...
    NoiseCorr_Analysis_dtheta(L23_Trials_evoked(1: Ne, :, :),...
    Pref_theta_E, dtheta_bound2, d_BinStep_micron, d_tot, 0);
fitpar_err_noi_E_dtheta2 = (CI_noi_E_dtheta2(:, 2) - CI_noi_E_dtheta2(:, 1)) / 2;
%
dtheta_bound_oppo = [60, 90] * (pi / 180);
[~, NoiseCorr_mean_E_dtheta_oppo, NoiseCorr_se_E_dtheta_oppo, fitpar_noi_E_dtheta_oppo, CI_noi_E_dtheta_oppo, ~] =...
    NoiseCorr_Analysis_dtheta(L23_Trials_evoked(1: Ne, :, :),...
    Pref_theta_E, dtheta_bound_oppo, d_BinStep_micron, d_tot, 0);
fitpar_err_noi_E_dtheta_oppo = (CI_noi_E_dtheta_oppo(:, 2) - CI_noi_E_dtheta_oppo(:, 1)) / 2;


% Raster example
load([pwd, '/Results/Results_Par', num2str(Parset_i), '_Simul1.mat'], 's_sample');


clear N_Trial_per_simul N_Simul N_trial_tot N_trial theta_stim



save([pwd, '/Results_Total_', num2str(Parset_i), '.mat']);

import java.lang.System
java.lang.System.exit(0)


