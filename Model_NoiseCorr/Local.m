%% Run for the first time
if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
addpath(genpath([dir_0, '/Study/CompNeuro/Projects/Functions_simul/']));
cd([dir_0, '/Study/CompNeuro/Projects/Micro-clustering']);
addpath(genpath([pwd, '/functions']));
addpath(genpath([pwd, '/Model_SigCorr/SpkSims']));    % Common functions
cd([pwd, '/Model_NoiseCorr']);
if ~exist([pwd, '/Parameters'], 'dir'), mkdir([pwd, '/Parameters']); end
if ~exist([pwd, '/Figures'], 'dir'), mkdir([pwd, '/Figures']); end
if ~exist([pwd, '/Results'], 'dir'), mkdir([pwd, '/Results']); end
% [pwd, '/functions']
addpath(genpath(pwd));
rng('shuffle');


% Parameters_FFWD;
% Use 10

% T_on = 400;

Parset_i = 1;
W = [3.2, 1.08, -3.5; 6.4, 0.8, -6.15];
K_in = repmat([1000 200 200], 2, 1);
sigma_micron = [144, 146, 110; 3.75, 20, 9.25];
kappa = [0.075, 0.125, 0.125];
alpha_slow = 0.3;    % !
Parameters_Recurrent_spatial(Parset_i, W, K_in, sigma_micron, kappa, alpha_slow);
%
Parset_i_dthetadep = 2; Parset_i_spatial = 1;
ifDthetaDep = [1, 1, 1; 1, 1, 1]; sigma_dtheta = 40 * (pi / 180); a_dtheta = 0.1;
Parameters_Recurrent_dthetadep(Parset_i_spatial, Parset_i_dthetadep, ifDthetaDep, sigma_dtheta, a_dtheta);


% kappa 0
Parset_i = 3;
W = [3.2, 1.08, -3.5; 6.4, 0.8, -6.15];
K_in = repmat([1000 200 200], 2, 1);
sigma_micron = [144, 146, 110; 11.4514, 11.4514, 11.4514];
kappa = [0, 0, 0];
alpha_slow = 0.3;
Parameters_Recurrent_spatial(Parset_i, W, K_in, sigma_micron, kappa, alpha_slow);


Parset_i_dthetadep = 5; Parset_i_spatial = 4;
ifDthetaDep = [1, 1, 1; 1, 1, 1]; sigma_dtheta = 40 * (pi / 180); a_dtheta = 0.1;
Parameters_Recurrent_dthetadep(Parset_i_spatial, Parset_i_dthetadep, ifDthetaDep, sigma_dtheta, a_dtheta);


% 4/6 uses Parameters_Recurrent_1/3.mat.
% 1 - 3 uses FFWD_spiking_generation.m, i.e. No real noise corr in L4
% 4 - 6 uses FFWD_spiking_generation_LIF.m, i.e. Noise corr in L4 due to J_matrix from X to F.
