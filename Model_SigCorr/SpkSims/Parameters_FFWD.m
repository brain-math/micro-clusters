if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
dir_L4 = [dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Model_SigCorr/L4_discrete_model'];


%% Exactly the same as main_build_up_L4 except dir name and some wrap-up
Nf_sqrt = 100; Nf = Nf_sqrt ^ 2;
d0_F = 7.5;    % interval distance between neighboring grids. in micron.
d_tot = Nf_sqrt * d0_F;

% Build up L4 tuning curves with not-bad tuning and corr. structure in the dataset.
% Define a imaginary layer "X", then filter it with a narrow Gaussian filter, to get L4 tuning curves.
Pref_theta_X_name = [dir_L4, '/Parameters/theta_gabor_SP_', num2str(Nf_sqrt), '.mat'];
if isfile(Pref_theta_X_name), load(Pref_theta_X_name, 'Pref_theta_X');
else, Pref_theta_X = pi * rand([Nf_sqrt, Nf_sqrt]); save(Pref_theta_X_name, 'Pref_theta_X'); end
%
N_stim = 10; theta_stim = [1: N_stim] * (pi / N_stim);
TN = @(r0, dtheta) r0 * exp(- (dtheta / (20 * (pi / 180))) .^ 2);    % Tuning curve shape in layer X.
%x = linspace(-pi/2, pi/2, 181); plot(x * (180/pi), TN(x)); xlim([-90 90]);
r_max = 10;
%
savename_conn_F = [dir_L4, '/Parameters/Conn_ffwd_', num2str(Nf_sqrt), '.mat'];
load(savename_conn_F, 'J_matrix');    % a sparse matrix
J_matrix_unnormalized = J_matrix; clear J_matrix

rX = NaN(Nf, N_stim);
rF_avg_unnormalized = NaN(Nf, N_stim);
for k = 1: N_stim
    rX_tmp = TN(r_max, dtheta(theta_stim(k), Pref_theta_X(:)));
    rX(:, k) = rX_tmp / 1000;    % kHz;
    rF_avg_unnormalized(:, k) = J_matrix_unnormalized * rX_tmp;
end
%rF_norm_factor = mean(max(rX * 1000, [], 2)) / mean(max(rF_avg_unnormalized, [], 2));
rF_norm_factor = r_max / mean(max(rF_avg_unnormalized, [], 2));
J_matrix = J_matrix_unnormalized * rF_norm_factor;
rF_avg = rF_avg_unnormalized * rF_norm_factor;    % Hz    % just for reference, won't be used in SpkSim.
% real rF(:, k) should be J_matrix * (rX(:, k) + noise);
clear k rX_tmp TN rF_norm_factor J_matrix_unnormalized rF_avg_unnormalized

% d_BinStep_micron = 5; N_sample = Nf;
% [gOSI_F, Pref_theta_F, d_BinCenter_F, SigCorr_mean_F, SigCorr_se_F,...
%     fitpar_exp2_F, CI_exp2_F, fitpar_exp1_F, CI_exp1_F] =...
%     FR_Analysis(rF_avg, theta_stim, d_BinStep_micron, d_tot, N_sample);

rF_off = 3e-3;    % kHz
%T_off = 100; T_on = 400;    % ms
sigma_n = 3 * 0.01;    % Input noise variance    % 3.5
tau_n = 40;    % Decay constant of input noise, ms.

clear dir_0 dir_L4
save([pwd, '/Parameters/Parameters_FFWD.mat']);

