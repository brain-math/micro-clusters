function [dist_ctr, JN_d_1d_pdf, JN_d_par, dtheta_ctr, JN_dtheta_pdf, JN_dtheta_par] =...
    Wiring_analysis(JN, N_post_sample, d_pre, Pref_theta_post, Pref_theta_pre)

%load([pwd, '/Parameters/JN_2.mat'], 'JN_ef', 'JN_ee', 'JN_ei', 'K_in',...
%    'kappa', 'sigma_micron', 'a_dtheta', 'sigma_dtheta');
%load([pwd, '/Parameters/Parameters_FFWD.mat'], 'Pref_theta_F');
%load([pwd, '/Results/Results_Total_2.mat'], 'Pref_theta_E', 'Pref_theta_I');
%
%load([pwd, '/Parameters/JN_2.mat'], 'JN_ef', 'K_in',...
%    'kappa', 'sigma_micron', 'a_dtheta', 'sigma_dtheta');
%load([pwd, '/Parameters/Parameters_FFWD.mat'], 'Pref_theta_F');
%load([pwd, '/Results/Results_Total_2.mat'], 'Pref_theta_E', 'Pref_theta_I');
% JN = JN_ef; clear JN_ef
% N_post_sample = 4000; d_pre = 7.5; Pref_theta_post = Pref_theta_E; Pref_theta_pre = Pref_theta_F;


N_post = size(JN, 1); N_pre = size(JN, 2);
N_sqrt_post = sqrt(N_post); N_sqrt_pre = sqrt(N_pre);
if N_post_sample > N_post, N_post_sample = N_post; end

d_circ = @(d, T) acos(cos(d * (2 * pi / T))) * (T / (2 * pi));
% figure; x = 0: 200; plot(x, acos(cos(2 * pi * x / 100)) * (50 / pi));
% axis([0 200 0 200]); axis square;
%
idx_post_sample = randperm(N_post); idx_post_sample = sort(idx_post_sample(1: N_post_sample));
K = sum(JN(1, :));
dist_sample = NaN(K * N_post_sample, 1);
JN_sample = NaN(K * N_post_sample, 1);
dtheta_sample = NaN(K * N_post_sample, 1);
icount = 0;
%
for k = 1: N_post_sample
    idx_post_k = idx_post_sample(k);
    idx_pre_connected = find(JN(idx_post_k, :)' > 0);
    N_idx_pre = length(idx_pre_connected);
    %
    i_post_on_pre = (mod(idx_post_k - 1, N_sqrt_post) + 1) * (N_sqrt_pre / N_sqrt_post);
    j_post_on_pre = ceil(idx_post_k / N_sqrt_post) * (N_sqrt_pre / N_sqrt_post);
    i_pre_connected = (mod(idx_pre_connected - 1, N_sqrt_pre) + 1);
    j_pre_connected = ceil(idx_pre_connected / N_sqrt_pre);
    dist_sample(icount + 1: icount + N_idx_pre) =...
        sqrt(d_circ(i_pre_connected - i_post_on_pre, N_sqrt_pre) .^ 2 +...
        d_circ(j_pre_connected - j_post_on_pre, N_sqrt_pre) .^ 2);
    %
    JN_sample(icount + 1: icount + N_idx_pre) = JN(idx_post_k, idx_pre_connected)';
    %
    dtheta_sample(icount + 1: icount + N_idx_pre) =...
        abs(dtheta(Pref_theta_post(idx_post_k), Pref_theta_pre(idx_pre_connected)));
    %
    icount = icount + N_idx_pre;
end
dist_sample = dist_sample(1: icount) * d_pre;
JN_sample = JN_sample(1: icount);
dtheta_sample = dtheta_sample(1: icount) * (180 / pi);    % deg


dist_BinStep = 5; dist_ctr_idx0 = 2;%3;
% Compare for multiple values: 2.5, 5, 7.5, 10, 15;
% Large BinStep makes far away curve smoother, but artifactually increases the scale of  narrow component.
% Finally best is 5, but start from 3rd point 12.5 -- which 2.5 and 5 overlaps, but failed for >= 7.5
dist_min = 0; dist_max = ceil(N_sqrt_pre * d_pre / (sqrt(2) * dist_BinStep)) * dist_BinStep;
dist_ctr = [dist_min + (dist_BinStep / 2): dist_BinStep: dist_max - (dist_BinStep / 2)]';
dist_edge = [dist_min: dist_BinStep: dist_max]'; 
%
[JN_d_avg, ~, JN_d_num] = histogram_mean_sem(JN_sample, dist_sample, dist_edge);
JN_d_1d = JN_d_avg .* JN_d_num;
JN_d_1d_pdf = JN_d_1d / (K * N_post_sample * dist_BinStep);
% sum(JN_d_1d_pdf) * dist_BinStep    % this integ. of pdf should be strictly 1
dist_ctr = dist_ctr(dist_ctr_idx0: end);
JN_d_1d_pdf = JN_d_1d_pdf(dist_ctr_idx0: end);
%JN_d_2d_pdf = JN_d_1d_pdf ./ (2 * pi * dist_ctr);
%plot(dist_ctr, JN_d_1d_pdf); hold on; xlim([0 300]);


dtheta_min = 0; dtheta_max = 90; dtheta_BinStep = 5;
dtheta_ctr = [dtheta_min + (dtheta_BinStep / 2): dtheta_BinStep: dtheta_max - (dtheta_BinStep / 2)]';
dtheta_edge = [dtheta_min: dtheta_BinStep: dtheta_max]';
%
[JN_dtheta_avg, ~, JN_dtheta_num] = histogram_mean_sem(JN_sample, dtheta_sample, dtheta_edge);
JN_dtheta_pdf_raw = JN_dtheta_avg .* JN_dtheta_num;
JN_dtheta_pdf = JN_dtheta_pdf_raw / (K * N_post_sample * dtheta_BinStep);
% sum(JN_dtheta_pdf) * dtheta_BinStep    % this integ. of pdf should be strictly 1




options = optimoptions('lsqnonlin', 'Display', 'none', 'MaxFunEvals', 1800,...
    'MaxIter', 1800, 'TolX', 1e-10, 'TolFun', 1e-10);
%
Gau2 = @(r, sigma) exp(- r .^ 2 / (2 * sigma ^ 2)) / (2 * pi * sigma ^ 2);
DoubleGau_2 = @(r, A, kappa, sigma_n, sigma_b)...
    A * ( kappa * Gau2(r, sigma_n) + (1 - kappa) * Gau2(r, sigma_b) );
DoubleGau_1 = @(r, A, kappa, sigma_n, sigma_b)...
    DoubleGau_2(r, A, kappa, sigma_n, sigma_b) .* (2 * pi * r);
DoubleGau_1_fit = @(r, A, kappa100, sigma_n10, sigma_b)...
    DoubleGau_1(r, A, kappa100 / 100, sigma_n10 / 10, sigma_b);
%
IniVal = [1, 0.1 * 100, 7.5 * 10, 135];
Err1 = @(par) DoubleGau_1_fit(dist_ctr, par(1), par(2), par(3), par(4)) - JN_d_1d_pdf;
JN_d_par = lsqnonlin(Err1, IniVal, zeros(1, 4), 1e10 * ones(1, 4), options);
JN_d_par(2) = JN_d_par(2) / 100; JN_d_par(3) = JN_d_par(3) / 10; 
% r0 = 0: 0.05: 500; plot(r0, DoubleGau_1(r0, JN_d_par(1),...
%     JN_d_par(2), JN_d_par(3), JN_d_par(4)), 'r');
% Overlap this with dist_BinStep = [2.5, 5, 7.5, 10, 15] !
%
%
Gau0 = @(dtheta, A_dtheta, sigma_dtheta, a_dtheta)...
    A_dtheta * exp(- dtheta.^2 / (2 * sigma_dtheta ^ 2)) + a_dtheta;
IniVal = [10e-3, 35, 6e-3];
Err2 = @(par) Gau0(dtheta_ctr(2: end), par(1), par(2), par(3)) - JN_dtheta_pdf(2: end);
JN_dtheta_par = lsqnonlin(Err2, IniVal, zeros(1, 3), 1e10 * ones(1, 3), options);
%
%figure; plot(dtheta_ctr, JN_dtheta_pdf); hold on; axis([0 90 0 0.02]);
%theta0 = 0: 0.1: 90; plot(theta0, Gau0(theta0, JN_dtheta_par(1), JN_dtheta_par(2), JN_dtheta_par(3)));









%N_post = size(JN_ef, 1);

%Pref_theta_post = Pref_theta_E(:);    % (N_post, 1)
%Pref_theta_pre = Pref_theta_F(:)';    % (1, N_pre)
%dtheta_post_pre = abs(bsxfun(@dtheta, Pref_theta_post, Pref_theta_pre));

%idx_conn = find(JN_ef > 0);

%dtheta_min = 0; dtheta_max = pi / 2; dtheta_BinStep = (dtheta_max - dtheta_min) / 100;
%dtheta_ctr = dtheta_min: dtheta_BinStep: dtheta_max;
%dtheta_ctr_deg = dtheta_ctr(2: end - 1) * (180 / pi);
%dtheta_edge = dtheta_min - (dtheta_BinStep / 2): dtheta_BinStep: dtheta_max + (dtheta_BinStep / 2);

%[JN_dtheta_avg, ~, JN_dtheta_num] = histogram_mean_sem(JN_ef(idx_conn), dtheta_post_pre(idx_conn), dtheta_edge);
%JN_dtheta_pdf_raw = JN_dtheta_avg .* JN_dtheta_num;
%JN_dtheta_pdf = JN_dtheta_pdf_raw(2: end - 1) / (K_in(1, 1) * N_post * dtheta_BinStep);

%figure; plot(dtheta_ctr_deg, JN_dtheta_pdf); axis([0 90 0 2]);

