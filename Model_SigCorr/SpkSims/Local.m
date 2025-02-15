%% Run for the first time
if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
addpath(genpath([dir_0, '/Study/CompNeuro/Projects/Functions_simul/']));
cd([dir_0, '/Study/CompNeuro/Projects/Micro-clustering']);
addpath(genpath([pwd, '/functions'])); cd([pwd, '/Model_SigCorr/SpkSims']);
if ~exist([pwd, '/Parameters'], 'dir'), mkdir([pwd, '/Parameters']); end
if ~exist([pwd, '/Figures'], 'dir'), mkdir([pwd, '/Figures']); end
if ~exist([pwd, '/Results'], 'dir'), mkdir([pwd, '/Results']); end
addpath(genpath(pwd));
rng('shuffle');


% Parameters_FFWD;


%% use 100, 7.5 FFWD, rF max 10.
% 1, 2 for d = 2.5 data
% 3, 4 for d = 5 data

Parset_i = 1;
W = [3.2, 1.08, -3.5; 6.4, 0.8, -6.15];
K_in = repmat([1000 200 200], 2, 1);
sigma_micron = [144, 146, 110; 3.75, 20, 9.25];
kappa = [0.075, 0.125, 0.125];
alpha_slow = 0.3;
Parameters_Recurrent_spatial(Parset_i, W, K_in, sigma_micron, kappa, alpha_slow);
%
Parset_i_dthetadep = 2; Parset_i_spatial = 1;
ifDthetaDep = [1, 1, 1; 1, 1, 1]; sigma_dtheta = 40 * (pi / 180); a_dtheta = 0.1;
Parameters_Recurrent_dthetadep(Parset_i_spatial, Parset_i_dthetadep, ifDthetaDep, sigma_dtheta, a_dtheta);



W = [3.2, 1.08, -3.5; 6.4, 0.8, -6.15];
K_in = repmat([1000 200 200], 2, 1);
%
% sigma_broad
Parset_i = 3;
sigma_micron = [100, 150, 150; 11.4514, 11.4514, 11.4514]; kappa = [0, 0, 0]; alpha_slow = 0;
Parameters_Recurrent_spatial(Parset_i, W, K_in, sigma_micron, kappa, alpha_slow);
%main_test(Parset_i);
%
Parset_i = 4;
sigma_micron = [150, 150, 150; 11.4514, 11.4514, 11.4514]; kappa = [0, 0, 0]; alpha_slow = 0;
Parameters_Recurrent_spatial(Parset_i, W, K_in, sigma_micron, kappa, alpha_slow);
%
Parset_i = 5;
sigma_micron = [150, 100, 100; 11.4514, 11.4514, 11.4514]; kappa = [0, 0, 0]; alpha_slow = 0;
Parameters_Recurrent_spatial(Parset_i, W, K_in, sigma_micron, kappa, alpha_slow);
%
% kappa
Parset_i = 6;
sigma_micron = [144, 146, 110; 4, 10, 10]; kappa = [0.025, 0.025, 0.025]; alpha_slow = 0;
Parameters_Recurrent_spatial(Parset_i, W, K_in, sigma_micron, kappa, alpha_slow);
%
Parset_i = 7;
sigma_micron = [144, 146, 110; 4, 10, 10]; kappa = [0.125, 0.025, 0.025]; alpha_slow = 0;
Parameters_Recurrent_spatial(Parset_i, W, K_in, sigma_micron, kappa, alpha_slow);
%
Parset_i = 8;
sigma_micron = [144, 146, 110; 4, 10, 10]; kappa = [0.125, 0.125, 0.125]; alpha_slow = 0;
Parameters_Recurrent_spatial(Parset_i, W, K_in, sigma_micron, kappa, alpha_slow);
%
% sigma_narrow
Parset_i = 9;
sigma_micron = [144, 146, 110; 10, 4, 4]; kappa = [0.075, 0.125, 0.125]; alpha_slow = 0;
Parameters_Recurrent_spatial(Parset_i, W, K_in, sigma_micron, kappa, alpha_slow);
%
Parset_i = 10;
sigma_micron = [144, 146, 110; 4, 4, 4]; kappa = [0.075, 0.125, 0.125]; alpha_slow = 0;
Parameters_Recurrent_spatial(Parset_i, W, K_in, sigma_micron, kappa, alpha_slow);
%
Parset_i = 11;
sigma_micron = [144, 146, 110; 4, 10, 10]; kappa = [0.075, 0.125, 0.125]; alpha_slow = 0;
Parameters_Recurrent_spatial(Parset_i, W, K_in, sigma_micron, kappa, alpha_slow);




W = [3.2, 1.08, -3.5; 6.4, 0.8, -6.15];
K_in = repmat([1000 200 200], 2, 1);
ifDthetaDep = [1, 1, 1; 1, 1, 1]; sigma_dtheta = 40 * (pi / 180); a_dtheta = 0.1;
%
for Parset_i_spatial = 3: 11
Parset_i_dthetadep = Parset_i_spatial + 9;
Parameters_Recurrent_dthetadep(Parset_i_spatial, Parset_i_dthetadep, ifDthetaDep, sigma_dtheta, a_dtheta);
end









