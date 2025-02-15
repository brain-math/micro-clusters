function getRecDthetaDepAdjacencyMatrices(savename_result_before, savename_JN_before,...
    ifDthetaDep, d_tot, sigma_dtheta, a_dtheta, savename_JN_after)
% Pref_theta_F, Pref_theta_E, Pref_theta_I must be in square
% sigma_dtheta should be in rad

load(savename_result_before, 'Pref_theta_F', 'Pref_theta_E', 'Pref_theta_I');

load(savename_JN_before, 'K_in', 'sigma_micron', 'kappa');

K_ef = K_in(1, 1); K_if = K_in(2, 1); kappa_f = kappa(1);
sigma_n_f = sigma_micron(2, 1); sigma_b_f = sigma_micron(1, 1);    % ef and if
if ifDthetaDep(1, 1) == 1, JN_ef =...
    getDthetaDepAdjacencyMatrix(d_tot, kappa_f, sigma_n_f, sigma_b_f,...
        K_ef, Pref_theta_E, Pref_theta_F, sigma_dtheta, a_dtheta);
else, load(savename_JN_before, 'JN_ef'); end
if ifDthetaDep(2, 1) == 1, JN_if =...
    getDthetaDepAdjacencyMatrix(d_tot, kappa_f, sigma_n_f, sigma_b_f,...
        K_if, Pref_theta_I, Pref_theta_F, sigma_dtheta, a_dtheta);
else, load(savename_JN_before, 'JN_if'); end
%JN = getDthetaDepAdjacencyMatrix(d_tot, kappa, sigma_n, sigma_b, K_in,...
%    Pref_theta_post, Pref_theta_pre, sigma_dtheta, a_dtheta)


K_ee = K_in(1, 2); K_ie = K_in(2, 2); kappa_e = kappa(2);
sigma_n_e = sigma_micron(2, 2); sigma_b_e = sigma_micron(1, 2);    % ee and ie
if ifDthetaDep(1, 2) == 1, JN_ee =...
    getDthetaDepAdjacencyMatrix(d_tot, kappa_e, sigma_n_e, sigma_b_e,...
        K_ee, Pref_theta_E, Pref_theta_E, sigma_dtheta, a_dtheta);
else, load(savename_JN_before, 'JN_ee'); end
if ifDthetaDep(2, 2) == 1, JN_ie =...
    getDthetaDepAdjacencyMatrix(d_tot, kappa_e, sigma_n_e, sigma_b_e,...
        K_ie, Pref_theta_I, Pref_theta_E, sigma_dtheta, a_dtheta);
else, load(savename_JN_before, 'JN_ie'); end


K_ei = K_in(1, 3); K_ii = K_in(2, 3); kappa_i = kappa(3);
sigma_n_i = sigma_micron(2, 3); sigma_b_i = sigma_micron(1, 3);    % ei and ii
if ifDthetaDep(1, 3) == 1, JN_ei =...
    getDthetaDepAdjacencyMatrix(d_tot, kappa_i, sigma_n_i, sigma_b_i,...
        K_ei, Pref_theta_E, Pref_theta_I, sigma_dtheta, a_dtheta);
else, load(savename_JN_before, 'JN_ei'); end
if ifDthetaDep(2, 3) == 1, JN_ii =...
    getDthetaDepAdjacencyMatrix(d_tot, kappa_i, sigma_n_i, sigma_b_i,...
        K_ii, Pref_theta_I, Pref_theta_I, sigma_dtheta, a_dtheta);
else, load(savename_JN_before, 'JN_ii'); end

save(savename_JN_after, 'JN_ef', 'JN_if', 'JN_ee', 'JN_ei', 'JN_ie', 'JN_ii',...
    'K_in', 'sigma_micron', 'kappa', 'ifDthetaDep', 'sigma_dtheta', 'a_dtheta');









