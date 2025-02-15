function JN = getDthetaDepAdjacencyMatrix(d_tot, kappa, sigma_n, sigma_b, K_in,...
    Pref_theta_post, Pref_theta_pre, sigma_dtheta, a_dtheta)
% Pref_theta_post, Pref_theta_pre must be in square
% sigma_dtheta should be in rad

N_sqrt_post = size(Pref_theta_post, 1); N_sqrt_pre = size(Pref_theta_pre, 1);
N_post = N_sqrt_post ^ 2; N_pre = N_sqrt_pre ^ 2;
d_pre = d_tot / N_sqrt_pre; d_post = d_tot / N_sqrt_post;

Coor_post_x = [d_post: d_post: d_tot] - (d_post / 2); Coor_post_y = Coor_post_x(end: -1: 1);
Coor_pre_x = [d_pre: d_pre: d_tot] - (d_pre / 2); Coor_pre_y = Coor_pre_x(end: -1: 1);

Gau2 = @(d, sigma) exp(- d .^ 2 / (2 * sigma ^ 2)) / (2 * pi * sigma ^ 2);
Gau2_mix = @(d, kappa, sigma_narrow, sigma_broad)...
    kappa * Gau2(d, sigma_narrow) + (1 - kappa) * Gau2(d, sigma_broad);
P_dtheta = @(dtheta, sigma, a) (1 - a) * exp(- dtheta .^ 2 / (2 * sigma ^ 2)) + a;
% d = -90: 0.1: 90; y = P_dtheta(d, 20, 0.18); figure; plot(d, y); axis([-90 90 0 1]);

rndset = rand(K_in, N_post);

JN = NaN(N_post, N_pre);

for k_post = 1: N_post
    i_post = mod(k_post - 1, N_sqrt_post) + 1; y_post = Coor_post_y(i_post);
    j_post = ceil(k_post / N_sqrt_post); x_post = Coor_post_x(j_post);
    delta_x = min(abs([Coor_pre_x - x_post; Coor_pre_x - (x_post - d_tot); Coor_pre_x - (x_post + d_tot)]), [], 1);
    delta_y = min(abs([Coor_pre_y - y_post; Coor_pre_y - (y_post - d_tot); Coor_pre_y - (y_post + d_tot)]), [], 1)';
    dist = sqrt(bsxfun(@plus, delta_x .^ 2, delta_y .^ 2));
    pdf_raw_dist = Gau2_mix(dist, kappa, sigma_n, sigma_b);
    dtheta_matrix = real(dtheta(Pref_theta_post(i_post, j_post), Pref_theta_pre));
    % I dont know fucking why + 0.0000i comes out
    pdf_raw_dtheta = P_dtheta(dtheta_matrix, sigma_dtheta, a_dtheta);
    pdf_raw = pdf_raw_dist .* pdf_raw_dtheta; pdf_raw = pdf_raw / max(pdf_raw(:));
    idx_pre = pdf_rand(pdf_raw, 0, rndset(:, k_post));
    JN(k_post, :) = histcounts(idx_pre, [0.5: 1: N_pre + 0.5]);
end

