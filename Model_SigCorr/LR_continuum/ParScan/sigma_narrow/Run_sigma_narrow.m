if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
addpath(genpath([dir_0, '/Study/CompNeuro/Projects/Functions_simul/']));
dir_local = [dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Model_SigCorr/LR_continuum'];
cd(dir_local); addpath([dir_local, '/Core']);
rng('shuffle');


load([dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging/',...
    'Analysis_5_JointStat_Overall/SigCorr_plot_d5.mat'], 'fitpar_exp2_data2', 'Exp2');
fitpar_exp2_L4 = fitpar_exp2_data2(2, :); clear fitpar_exp2_data2    % L4 input, use exp2 model
%
d_micron = 0: 0.1: 100;
W = [3.2, 1.08, -3.5; 6.4, 0.8, -6.15];
sigma_micron = [144, 146, 110; NaN, NaN, NaN]; sigma_micron_B = sigma_micron(1, :);
kappa = [0.075, 0.125, 0.125];
%
scan_mode = 1;    % 1: e = i v.s. f; 2: i v.s. f = e
if scan_mode == 1
    y_name = '\sigma\_Narrow\_E = \sigma\_Narrow\_I (\mum)'; x_name = '\sigma\_Narrow\_F (\mum)';
elseif scan_mode == 2
    y_name = '\sigma\_Narrow\_I (\mum)'; x_name = '\sigma\_Narrow\_F = \sigma\_Narrow\_E (\mum)';
end
y_list = 2: 1: 15; Ny = length(y_list);
x_list = 2: 1: 15; Nx = length(x_list);
%
Cov_L23E = NaN(Ny, Nx, length(d_micron));
%Exp2 = @(x, A, lambda, b) A * exp(- (x / lambda) .^ 2) + b;
par_lb = [0, 0, 0]; par_ub = [1e5, 1e5, 1e5];
options = optimoptions('lsqnonlin', 'Display', 'none', 'MaxFunEvals', 1200,...
    'MaxIter', 1200, 'TolX', 1e-8, 'TolFun', 1e-8);
A_amp = 1e3;
Cov_L23E_par = NaN(Ny, Nx, 3);


for i = 1: Ny
for j = 1: Nx
    if scan_mode == 1
        sigma_micron(2, 2: 3) = y_list(i);    % sigma_N_E = sigma_N_I
        sigma_micron(2, 1) = x_list(j);    % sigma_N_F
    elseif scan_mode == 2
        sigma_micron(2, 3) = y_list(i);    % sigma_N_I
        sigma_micron(2, 1: 2) = x_list(j);    % sigma_N_F = sigma_N_E
    end
    [~, ~, Cov_L23E_ij, ~] = LR_continuum_FT_func(d_micron, fitpar_exp2_L4, W, sigma_micron, kappa);
    Cov_L23E(i, j, :) = reshape(Cov_L23E_ij, [1 1 length(d_micron)]);
    %
    Err = @(par) Exp2(d_micron, par(1), par(2), par(3)) - A_amp * Cov_L23E_ij;
    IniVal = [A_amp * max(Cov_L23E_ij), 10, 0];
    try
        par_tmp = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
        par_tmp([1 3]) = par_tmp([1 3]) / A_amp; par_tmp = reshape(par_tmp, [1 1 3]);
        Cov_L23E_par(i, j, :) = par_tmp;
    catch
        fprintf(['sigma_narrow_x = ', num2str(x_list(j)), ', sigma_narrow_y = ', num2str(y_list(i)), ' failed! \n']);
    end
end
fprintf([num2str(i), ' / ', num2str(Ny), ' end.\n']);
end
clear dir_0 sigma_micron par_lb par_ub options A_amp i j Cov_L23E_ij Err IniVal par_tmp 
%
save([dir_local, '/ParScan/sigma_narrow/sigma_narrow_mode', num2str(scan_mode), '.mat'],...
    'Cov_L23E', 'Cov_L23E_par', 'd_micron', 'Exp2', 'fitpar_exp2_L4', 'kappa',...
    'Nx', 'Ny', 'sigma_micron_B', 'W', 'x_list', 'x_name', 'y_list', 'y_name');

