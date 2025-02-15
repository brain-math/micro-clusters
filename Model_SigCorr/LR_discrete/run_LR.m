function run_LR(J0, Ge, Gi, savename_JN, savename_L4, ifIsyn, ifBroadOnly, d_BinStep_E_micron, savename_result)

load(savename_JN);
% 'JN_ef', 'JN_if', 'JN_ee', 'JN_ei', 'JN_ie', 'JN_ii', 'K_in', 'sigma_micron', 'kappa'
% 'ifDthetaDep', 'sigma_dtheta', 'a_dtheta'
Ne = size(JN_ef, 1); Ni = size(JN_if, 1); Nrec = Ne + Ni;
%
%J0 = W ./ (repmat([Ge; Gi], 1, 3) .* K_in);
% PSP (post synatic potential) of individual synapse, (V)
% [ef, ee, ei; if, ie, ii]
%
W_r = [Ge * J0(1, 2) * JN_ee, Ge * J0(1, 3) * JN_ei;...
    Gi * J0(2, 2) * JN_ie, Gi * J0(2, 3) * JN_ii];
clear JN_ee JN_ei JN_ie JN_ii
W_h = [Ge * J0(1, 1) * JN_ef; Gi * J0(2, 1) * JN_if];
clear JN_ef JN_if

load(savename_L4, 'N_stim', 'rF');
rR = NaN(Nrec, N_stim);
if ifIsyn == 1
    IsynF = NaN(Nrec, N_stim); IsynE = NaN(Nrec, N_stim); IsynI = NaN(Nrec, N_stim);    % positive IsynI
end
for stim_k = 1: N_stim
    h_k = W_h * rF(:, stim_k);
    rR_tmp = (eye(Nrec) - W_r) \ h_k;
    if ifIsyn == 1
        IsynF(:, Feature_i) = h_k;
        IsynE(:, Feature_i) = W_r(:, 1: Ne) * rR_tmp(1: Ne);
        IsynI(:, Feature_i) = - W_r(:, Ne + 1: Nrec) * rR_tmp(Ne + 1: Nrec);
    end
    rR(:, stim_k) = rR_tmp .* (rR_tmp > 0);    % ReLU
    fprintf([num2str(stim_k), ' / ', num2str(N_stim), '.\n']);
end
clear stim_k h_k rR_tmp W_r W_h

load(savename_L4, 'theta_stim', 'd_tot', 'Pref_theta_F');
if ifBroadOnly == 0

    [gOSI_E, Pref_theta_E, d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E,...
        fitpar_exp2_E, CI_exp2_E, fitpar_exp1_E, CI_exp1_E] =...
        FR_Analysis(rR(1: Ne, :), theta_stim, 0, d_tot, Ne);
    [gOSI_I, Pref_theta_I, d_BinCenter_I, SigCorr_mean_I, SigCorr_se_I,...
        fitpar_exp2_I, CI_exp2_I, fitpar_exp1_I, CI_exp1_I] =...
        FR_Analysis(rR(Ne + 1: end, :), theta_stim, 0, d_tot, Ni);
elseif ifBroadOnly == 1
    d_BinStep_E_micron = 15; d_BinStep_I_micron = 15;
    [gOSI_E, Pref_theta_E, d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E] =...
        FR_Analysis_NoFitting(rR(1: Ne, :), theta_stim, d_BinStep_E_micron, d_tot, Ne);
    [gOSI_I, Pref_theta_I, d_BinCenter_I, SigCorr_mean_I, SigCorr_se_I] =...
        FR_Analysis_NoFitting(rR(Ne + 1: end, :), theta_stim, d_BinStep_I_micron, d_tot, Ni);
end

variableList = {'rR', 'IsynF', 'IsynE', 'IsynI',...
    'sigma_micron', 'kappa', 'ifDthetaDep', 'sigma_dtheta', 'a_dtheta',...
    'K_in', 'Ge', 'Gi', 'J0', 'Pref_theta_F',...
    'gOSI_E', 'Pref_theta_E', 'd_BinCenter_E', 'SigCorr_mean_E', 'SigCorr_se_E',...
    'fitpar_exp2_E', 'CI_exp2_E', 'fitpar_exp1_E', 'CI_exp1_E',...
    'gOSI_I', 'Pref_theta_I', 'd_BinCenter_I', 'SigCorr_mean_I', 'SigCorr_se_I',...
    'fitpar_exp2_I', 'CI_exp2_I', 'fitpar_exp1_I', 'CI_exp1_I'};
for variableIndex = 1: length(variableList)
if exist(variableList{variableIndex}, 'var')
    if ~exist(savename_result, 'file')
        save(savename_result, variableList{variableIndex});
    else
        save(savename_result, variableList{variableIndex}, '-append');
    end
end
end
