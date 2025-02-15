function Results_Processing

Parset_i = str2num(getenv('SLURM_ARRAY_TASK_ID'));

load([pwd, '/Parameters/Parameters_FFWD.mat'], 'Nf', 'Nf_sqrt', 'd_tot', 'N_stim', 'theta_stim');
load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param');
Ne = param.Ne; Ne_sqrt = sqrt(Ne);
Ni = param.Ni; Ni_sqrt = sqrt(Ni); Nrec = param.N;
N_SpkCountperSimul = param.N_SpkCountperSimul; clear param
load([pwd, '/Results/Results_Par', num2str(Parset_i), '_Feature1_Repeat1.mat'], 'N_trial');    % N_trial = 200;
N_repeat = N_trial / N_SpkCountperSimul; clear N_SpkCountperSimul
%
L4_FR_tot = zeros(Nf, N_stim); L23_FR_tot = zeros(Nrec, N_stim); L23_Isyn = zeros(Nrec, N_stim, 4);    % F, E, I, Tot 
load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param');
K_in = param.K_in; sigma_micron = param.sigma_micron; kappa = param.kappa;
if isfield(param, 'ifDthetaDep'), ifDthetaDep = param.ifDthetaDep; end
if isfield(param, 'sigma_dtheta'), sigma_dtheta = param.sigma_dtheta; end
if isfield(param, 'a_dtheta'), a_dtheta = param.a_dtheta; end
alpha_slow = param.alpha_slow;
clear param
%
for Feature_i = 1: N_stim
    L4_FR_tmp = zeros(Nf, 1); L23_FR_tmp = zeros(Nrec, 1); Isyn_tot_tmp = zeros(Nrec, 1, 4); 
    for Repeat_i = 1: N_repeat
        load([pwd, '/Results/Results_Par', num2str(Parset_i), '_Feature',...
            num2str(Feature_i), '_Repeat', num2str(Repeat_i), '.mat'], 'L4_FR', 'L23_FR', 'Isyn');
        % (Nrec, N_stim)    % (Nrec, N_stim, 3)
        L4_FR_tmp = L4_FR_tmp + sum(L4_FR, 2) / N_trial;
        L23_FR_tmp = L23_FR_tmp + sum(L23_FR, 2) / N_trial;
        Isyn = cat(3, Isyn, sum(Isyn, 3));    % (Nrec, N_stim, 4)
        Isyn_tot_tmp = Isyn_tot_tmp + sum(Isyn, 2) / N_trial;
    end
    L4_FR_tot(:, Feature_i) = L4_FR_tmp;
    L23_FR_tot(:, Feature_i) = L23_FR_tmp;
    L23_Isyn(:, Feature_i, :) = Isyn_tot_tmp;
end
clear N_trial N_repeat Feature_i Repeat_i L4_FR_tmp L23_FR_tmp Isyn_tot_tmp L4_FR L23_FR Isyn
%
L4_FR = L4_FR_tot; clear L4_FR_tot
L23_FR = L23_FR_tot; clear L23_FR_tot
%
N_sample = [Nf, 5000, 2000];
d_BinStep_F_micron = 5;
[gOSI_F, Pref_theta_F, d_BinCenter_F, SigCorr_mean_F, SigCorr_se_F,...
    fitpar_exp2_F, CI_exp2_F, fitpar_exp1_F, CI_exp1_F] =...
    FR_Analysis(L4_FR, theta_stim, d_BinStep_F_micron, d_tot, N_sample(1));
%
if any(Parset_i == [3 4 5 12 13 14])    % Broad only
    d_BinStep_E_micron = 15; d_BinStep_I_micron = 15; 
    [gOSI_E, Pref_theta_E, d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E] = ...
        FR_Analysis_NoFitting(L23_FR(1: Ne, :), theta_stim, d_BinStep_E_micron, d_tot, N_sample(2));
    [gOSI_I, Pref_theta_I, d_BinCenter_I, SigCorr_mean_I, SigCorr_se_I] = ...
        FR_Analysis_NoFitting(L23_FR(Ne + 1: Nrec, :), theta_stim, d_BinStep_I_micron, d_tot, N_sample(3));
else
    d_BinStep_E_micron = 0; d_BinStep_I_micron = 0;    % Real coordinate before 51 and 10 after that.
    [gOSI_E, Pref_theta_E, d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E,...
        fitpar_exp2_E, CI_exp2_E, fitpar_exp1_E, CI_exp1_E] =...
        FR_Analysis(L23_FR(1: Ne, :), theta_stim, d_BinStep_E_micron, d_tot, N_sample(2));
    [gOSI_I, Pref_theta_I, d_BinCenter_I, SigCorr_mean_I, SigCorr_se_I] = ...
        FR_Analysis_NoFitting(L23_FR(Ne + 1: Nrec, :), theta_stim, d_BinStep_I_micron, d_tot, N_sample(3));
end
%
load([pwd, '/Results/Results_Par', num2str(Parset_i), '_Feature',...
    num2str(N_stim), '_Repeat1.mat'], 's_sample');


save([pwd, '/Results_Total_', num2str(Parset_i), '.mat']);

import java.lang.System
java.lang.System.exit(0)

