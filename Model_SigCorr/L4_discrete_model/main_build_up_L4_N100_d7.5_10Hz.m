clc; clear;
is_continuous = 0;    % Use a matrix from X to F (sampling from Gau. fiter), rather than imgaussfilt.

if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
addpath(genpath([dir_0, '/Study/CompNeuro/Projects/Functions_simul/']));
cd([dir_0, '/Study/CompNeuro/Projects/Micro-clustering']);
addpath(genpath([pwd, '/functions'])); 
cd([pwd, '/Model_SigCorr/L4_discrete_model']);
if ~exist([pwd, '/Parameters'], 'dir'), mkdir([pwd, '/Parameters']); end
addpath(genpath(pwd));
rng('shuffle');

Nf_sqrt = 100; Nf = Nf_sqrt ^ 2;
d0_F = 7.5;    % interval distance between neighboring grids. in micron.
d_tot = Nf_sqrt * d0_F;

% Build up L4 tuning curves with not-bad tuning and corr. structure in the dataset.
% Define a imaginary layer "X", then filter it with a narrow Gaussian filter, to get L4 tuning curves.
Pref_theta_X_name = [pwd, '/Parameters/theta_gabor_SP_', num2str(Nf_sqrt), '.mat'];
if isfile(Pref_theta_X_name), load(Pref_theta_X_name, 'Pref_theta_X');
else, Pref_theta_X = pi * rand([Nf_sqrt, Nf_sqrt]); save(Pref_theta_X_name, 'Pref_theta_X'); end
%
N_stim = 10; theta_stim = [1: N_stim] * (pi / N_stim);
TN = @(r0, dtheta) r0 * exp(- (dtheta / (20 * (pi / 180))) .^ 2);    % Tuning curve shape in layer X.
%x = linspace(-pi/2, pi/2, 181); plot(x * (180/pi), TN(x)); xlim([-90 90]);
r_max = 10;
%
sigma_F = 5.15;%5.8;
if is_continuous == 0
    K_F = 60;%25;
    pdf_Thr_F = 1e-4;    % use 0 will lead to repeat in cdf and error.
    savename_conn_F = [pwd, '/Parameters/Conn_ffwd_', num2str(Nf_sqrt), '.mat'];
    getConnGauF(Nf_sqrt, K_F, sigma_F, d_tot, pdf_Thr_F, savename_conn_F);
    load(savename_conn_F, 'J_matrix');
end

rX = NaN(Nf_sqrt, Nf_sqrt, N_stim);
rF = NaN(Nf_sqrt, Nf_sqrt, N_stim);
for k = 1: N_stim
    rX_tmp = TN(r_max, dtheta(theta_stim(k), Pref_theta_X));
    rX(:, :, k) = rX_tmp;
    if is_continuous == 1, rF(:, :, k) = imgaussfilt(rX_tmp, sigma_F / d0_F, 'Padding', 'circular');
    elseif is_continuous == 0, rF(:, :, k) = reshape(J_matrix * reshape(rX_tmp, [Nf, 1]), [Nf_sqrt, Nf_sqrt]); end
end
r_max_F = max(rF, [], 3); r_max_F = mean(r_max_F(:));
rF = reshape(rF * (r_max / r_max_F), [Nf, N_stim]);
% [mean(rF(:)), mean(max(rF, [], 2))]
clear k rX_tmp r_max_X r_max_F J_matrix

% Check L4 tuning curves
% figure;
% idx = randi(Nf); plot([0 theta_stim * (180 / pi)], [rF(idx, end), rF(idx, :)]); axis([0 180 0 15]);
% Check spatial clustering of firing rates
% figure; imagesc(reshape(rF(:, randi(10)), [Nf_sqrt, Nf_sqrt])); colormap('gray'); colorbar; axis square; caxis([0 11]);


% Check L4 signal corr.
% This is for LR. For SpkSim, only the X -> F matrix and Pref_theta_X will be used and Stats of Ffwd will be re-calculated.
d_BinStep_micron = 0; N_sample = Nf;
[gOSI_F, Pref_theta_F, d_BinCenter_F, SigCorr_mean_F, SigCorr_se_F,...
    fitpar_exp2_F, CI_exp2_F, fitpar_exp1_F, CI_exp1_F] =...
    FR_Analysis(rF, theta_stim, d_BinStep_micron, d_tot, N_sample);

% Plot test
Exp_type = 2;
plotL4(d_BinCenter_F, SigCorr_mean_F, SigCorr_se_F,...
    fitpar_exp2_F, CI_exp2_F, fitpar_exp1_F, CI_exp1_F, Exp_type);

%
pause(2); print(gcf, '-dpng', [pwd, '/plotL4_', num2str(Nf_sqrt), '.png']);
close;


if is_continuous == 0
save([pwd, '/Results_L4_', num2str(Nf_sqrt), '.mat'],...
    'Nf_sqrt', 'Nf', 'd0_F', 'd_tot', 'N_stim', 'theta_stim', 'TN', 'r_max',...
    'sigma_F', 'K_F', 'pdf_Thr_F', 'savename_conn_F', 'rX', 'rF', 'd_BinStep_micron', 'N_sample',...
    'gOSI_F', 'Pref_theta_F', 'd_BinCenter_F', 'SigCorr_mean_F', 'SigCorr_se_F',...
    'fitpar_exp2_F', 'CI_exp2_F', 'fitpar_exp1_F', 'CI_exp1_F');
elseif is_continuous == 1
save([pwd, '/Results_L4_', num2str(Nf_sqrt), '.mat'],...
    'Nf_sqrt', 'Nf', 'd0_F', 'd_tot', 'N_stim', 'theta_stim', 'TN', 'r_max',...
    'sigma_F', 'rX', 'rF', 'd_BinStep_micron', 'N_sample',...
    'gOSI_F', 'Pref_theta_F', 'd_BinCenter_F', 'SigCorr_mean_F', 'SigCorr_se_F',...
    'fitpar_exp2_F', 'CI_exp2_F', 'fitpar_exp1_F', 'CI_exp1_F');
end


