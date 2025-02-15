function getRecGauMixAdjacencyMatrices(d_tot, Nf_sqrt, Ne_sqrt, Ni_sqrt, K_in, sigma_micron, kappa, savename_JN)
% Example:
% sigma_micron = [144, 146, 110; 3.5, 22.5, 8.75]; kappa = [0.05, 0.05, 0.125];
% K_in = repmat([1000, 1000, 750], 2, 1);

d_f = d_tot / Nf_sqrt; d_e = d_tot / Ne_sqrt; d_i = d_tot / Ni_sqrt;

K_ef = K_in(1, 1); K_if = K_in(2, 1); kappa_f = kappa(1);
sigma_n_f = sigma_micron(2, 1); sigma_b_f = sigma_micron(1, 1);    % ef and if
JN_ef = getGauMixAdjacencyMatrix(Ne_sqrt, Nf_sqrt, d_f, kappa_f, sigma_n_f, sigma_b_f, K_ef);
JN_if = getGauMixAdjacencyMatrix(Ni_sqrt, Nf_sqrt, d_f, kappa_f, sigma_n_f, sigma_b_f, K_if);
% (N_sqrt_post, N_sqrt_pre, d_pre, kappa, sigma_n, sigma_b, K_in)

K_ee = K_in(1, 2); K_ie = K_in(2, 2); kappa_e = kappa(2);
sigma_n_e = sigma_micron(2, 2); sigma_b_e = sigma_micron(1, 2);    % ee and ie
JN_ee = getGauMixAdjacencyMatrix(Ne_sqrt, Ne_sqrt, d_e, kappa_e, sigma_n_e, sigma_b_e, K_ee);
JN_ie = getGauMixAdjacencyMatrix(Ni_sqrt, Ne_sqrt, d_e, kappa_e, sigma_n_e, sigma_b_e, K_ie);
% (N_sqrt_post, N_sqrt_pre, d_pre, kappa, sigma_n, sigma_b, K_in)

K_ei = K_in(1, 3); K_ii = K_in(2, 3); kappa_i = kappa(3);
sigma_n_i = sigma_micron(2, 3); sigma_b_i = sigma_micron(1, 3);    % ei and ii
JN_ei = getGauMixAdjacencyMatrix(Ne_sqrt, Ni_sqrt, d_i, kappa_i, sigma_n_i, sigma_b_i, K_ei);
JN_ii = getGauMixAdjacencyMatrix(Ni_sqrt, Ni_sqrt, d_i, kappa_i, sigma_n_i, sigma_b_i, K_ii);
% (N_sqrt_post, N_sqrt_pre, d_pre, kappa, sigma_n, sigma_b, K_in)

save(savename_JN, 'JN_ef', 'JN_if', 'JN_ee', 'JN_ei', 'JN_ie', 'JN_ii', 'K_in', 'sigma_micron', 'kappa');
