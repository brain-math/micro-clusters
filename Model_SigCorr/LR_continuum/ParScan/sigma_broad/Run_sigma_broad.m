if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
addpath(genpath([dir_0, '/Study/CompNeuro/Projects/Functions_simul/']));
dir_local = [dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Model_SigCorr/LR_continuum'];
cd(dir_local); addpath([dir_local, '/Core']);
rng('shuffle');


load([dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging/',...
    'Analysis_5_JointStat_Overall/SigCorr_plot_d5.mat'], 'fitpar_exp2_data2', 'Exp2');
fitpar_exp2_L4 = fitpar_exp2_data2(2, :); clear fitpar_exp2_data2    % L4 input, use exp2 model
%
d_micron = 0: 1: 500;
n_cycmicron = linspace(0, 0.2, 5001);
W = [3.2, 1.08, -3.5; 6.4, 0.8, -6.15];
kappa = [0, 0, 0];
%
scan_mode = 1;    % 1: e = i v.s. f; 2: i v.s. f = e
if scan_mode == 1
    y_name = '\sigma\_Broad\_E = \sigma\_Broad\_I (\mum)'; x_name = '\sigma\_Broad\_F (\mum)';
elseif scan_mode == 2
    y_name = '\sigma\_Broad\_I (\mum)'; x_name = '\sigma\_Broad\_F = \sigma\_Broad\_E (\mum)';
end
y_list = 75: 5: 175; Ny = length(y_list);
x_list = 75: 5: 175; Nx = length(x_list);
%
Cov_L23E = NaN(Ny, Nx, length(d_micron));
Cov_L23E_par = NaN(Ny, Nx, 2);
% 1: feature wavelength (\mum), 1 / n_max; 2: stength, peak value - 0 value
% If peak is at n = 0, then remain NaN


for i = 1: Ny
for j = 1: Nx
    sigma_micron = ones(2, 3);
    if scan_mode == 1
        sigma_micron(1, 2: 3) = y_list(i);    % sigma_B_E = sigma_B_I
        sigma_micron(1, 1) = x_list(j);    % sigma_B_F
    elseif scan_mode == 2
        sigma_micron(1, 3) = y_list(i);    % sigma_B_I
        sigma_micron(1, 1: 2) = x_list(j);    % sigma_B_F = sigma_B_E
    end
    [FR_n_L23E_ij, ~, Cov_L23E_ij, ~] = LR_continuum_FT_func(d_micron, fitpar_exp2_L4, W, sigma_micron, kappa);
    Cov_L23E(i, j, :) = reshape(Cov_L23E_ij, [1 1 length(d_micron)]);
    %
    [r_max, n_max] = max(FR_n_L23E_ij);
    if n_max ~= 1
        Cov_L23E_par(i, j, 1) = 1 / n_cycmicron(n_max);
        Cov_L23E_par(i, j, 2) = r_max - FR_n_L23E_ij(1);
    end
end
fprintf([num2str(i), ' / ', num2str(Ny), ' end.\n']);
end
clear dir_0 sigma_micron i j FR_n_L23E_ij Cov_L23E_ij r_max n_max
%
save([dir_local, '/ParScan/sigma_broad/sigma_broad_mode', num2str(scan_mode), '.mat'],...
    'Cov_L23E', 'Cov_L23E_par', 'd_micron', 'Exp2', 'fitpar_exp2_L4', 'kappa',...
    'n_cycmicron', 'Nx', 'Ny', 'W', 'x_list', 'x_name', 'y_list', 'y_name');


