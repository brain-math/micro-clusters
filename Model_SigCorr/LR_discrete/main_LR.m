%% Run for the first time
if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
addpath(genpath([dir_0, '/Study/CompNeuro/Projects/Functions_simul/']));
dir_local = [dir_0, '/Study/CompNeuro/Projects/Micro-clustering'];
cd(dir_local); addpath(genpath([dir_local, '/functions'])); cd([dir_local, '/Model_SigCorr/LR_discrete']);
if ~exist([pwd, '/Parameters'], 'dir'), mkdir([pwd, '/Parameters']); end
if ~exist([pwd, '/Figures'], 'dir'), mkdir([pwd, '/Figures']); end
if ~exist([pwd, '/Results'], 'dir'), mkdir([pwd, '/Results']); end
addpath(genpath(pwd));
rng('shuffle');


%% COMMON PARAMETERS (Run for all modules)
Nf_sqrt = 100; d0_F = 7.5;    % rF max = 10.
savename_L4 = [dir_local, '/Model_SigCorr/L4_discrete_model/Results_L4_', num2str(Nf_sqrt), '.mat'];
savename_JN_func = @(module_name) [pwd, '/Parameters/JN_', module_name, '.mat'];
savename_result_func = @(module_name) [pwd, '/Results/Results_', module_name, '.mat'];
savename_J0_func = @(module_name) [pwd, '/Parameters/J0_', module_name, '.mat'];
savename_figure_func = @(module_name) [pwd, '/Figures/', module_name];
%
Ne_sqrt = Nf_sqrt; Ni_sqrt = Ne_sqrt / 2;
load(savename_L4, 'd_tot');
ifIsyn = 0;
d_BinStep_E_micron = 5;    % um
%
K_in = repmat([1000 200 200], 2, 1);
W_est = [3.2, 1.08, -3.5; 6.4, 0.8, -6.15];
% Wii Must > 6, but 6.3 make it easier to explode.
% |Wii| > |Wei| * (Wif / Wef) - 1 = 3.5 * 2 - 1 = 6
% (eye(2) - W(:, 2: 3)) \ W(:, 1)    % [0.22; 0.92]
G_est = [12 * ones(1, 3); 15 * ones(1, 3)];    % V^(-1)
% K_in_ij * J0_ij = Jtot_ij = W_ij / G_ij
J0 = (W_est ./ G_est) ./ K_in;    % (V)




%% Pure spatial model
%module_name = 'Spatial_test_0'; ifBroadOnly = 0;
%
%K_in = repmat([600 600 400], 2, 1);
%%K_in_ij * J0_ij = Jtot_ij = W_ij / G_ij
%J0 = (W_est ./ G_est) ./ K_in;    % (V)
%Ge = 12; Gi = 15;
%sigma_micron = [144, 146, 110; 3.75, 15, 9.25];
%kappa = [0.075, 0.125, 0.125];
%
%getRecGauMixAdjacencyMatrices(d_tot, Nf_sqrt, Ne_sqrt, Ni_sqrt,...
%    K_in, sigma_micron, kappa, savename_JN_func(module_name));
%run_LR(J0, Ge, Gi, savename_JN_func(module_name), savename_L4,...
%    ifIsyn, ifBroadOnly, d_BinStep_E_micron, savename_result_func(module_name));
%plot_spatial(savename_L4, savename_result_func(module_name),...
%    ifBroadOnly, 2, savename_figure_func(module_name));




%% Dtheta dep
Ge_spatial = 12; Gi_spatial = 15;    % Gain_3.png    % 7, 27.5
Ge_ftrspc = 12; Gi_ftrspc = 15;    % Gain_4.png    % 10.7, 27.5
%
ifDthetaDep = [1, 1, 1; 1, 1, 1]; sigma_dtheta = 40 * (pi / 180); a_dtheta = 0.1;



module_name = 'DthetaDep_test_0'; ifBroadOnly = 0;
sigma_micron = [144, 146, 110; 3.75, 20, 9.25];
kappa = [0.075, 0.125, 0.125];
%
%% kappa
%module_name = 'DthetaDep_test_1'; ifBroadOnly = 0;
%sigma_micron = [144, 146, 110; 4, 10, 10];
%kappa = [0.025, 0.025, 0.025];
%%
%module_name = 'DthetaDep_test_2'; ifBroadOnly = 0;
%sigma_micron = [144, 146, 110; 4, 10, 10];
%kappa = [0.125, 0.025, 0.025];
%%
%module_name = 'DthetaDep_test_3'; ifBroadOnly = 0;
%sigma_micron = [144, 146, 110; 4, 10, 10];
%kappa = [0.125, 0.125, 0.125];
%%
%% sigma_narrow
%module_name = 'DthetaDep_test_4'; ifBroadOnly = 0;
%sigma_micron = [144, 146, 110; 10, 4, 4];
%kappa = [0.075, 0.125, 0.125];
%%
%module_name = 'DthetaDep_test_5'; ifBroadOnly = 0;
%sigma_micron = [144, 146, 110; 4, 4, 4];
%kappa = [0.075, 0.125, 0.125];
%%
%module_name = 'DthetaDep_test_6'; ifBroadOnly = 0;
%sigma_micron = [144, 146, 110; 4, 10, 10];
%kappa = [0.075, 0.125, 0.125];
%
%
getRecGauMixAdjacencyMatrices(d_tot, Nf_sqrt, Ne_sqrt, Ni_sqrt,...
    K_in, sigma_micron, kappa, savename_JN_func([module_name, '_before']));
run_LR(J0, Ge_spatial, Gi_spatial, savename_JN_func([module_name, '_before']), savename_L4,...
    ifIsyn, ifBroadOnly, d_BinStep_E_micron, savename_result_func([module_name, '_before']));
%
getRecDthetaDepAdjacencyMatrices(savename_result_func([module_name, '_before']),...
    savename_JN_func([module_name, '_before']), ifDthetaDep,...
    d_tot, sigma_dtheta, a_dtheta, savename_JN_func([module_name, '_after']));
run_LR(J0, Ge_ftrspc, Gi_ftrspc, savename_JN_func([module_name, '_after']), savename_L4,...
    ifIsyn, ifBroadOnly, d_BinStep_E_micron, savename_result_func([module_name, '_after']));
plot_dthetadep(savename_L4, savename_result_func([module_name, '_before']),...
    savename_result_func([module_name, '_after']), savename_JN_func([module_name, '_after']),...
    ifBroadOnly, 2, savename_figure_func(module_name));


