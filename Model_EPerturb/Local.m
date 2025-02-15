%% Run for the first time
if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
addpath(genpath([dir_0, '/Study/CompNeuro/Projects/Functions_simul']));
dir_local = [dir_0, '/Study/CompNeuro/Projects/Micro-clustering'];
cd(dir_local); addpath(genpath([pwd, '/functions'])); cd([pwd, '/Model_EPerturb']);
if ~exist([pwd, '/Parameters'], 'dir'), mkdir([pwd, '/Parameters']); end
if ~exist([pwd, '/Figures'], 'dir'), mkdir([pwd, '/Figures']); end
if ~exist([pwd, '/Results'], 'dir'), mkdir([pwd, '/Results']); end
addpath(genpath(pwd));
rng('shuffle');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 15000 trials of 2s off + 1s on (use last 500ms for FR), 10Hz stim
% 30% slow PSP; rX = 7.5Hz;
% sigma_micron = [144, 146, 110; 3.75, 15 or 20, 9.25]; kappa = [0.075, 0.125, 0.125];
% K_in = [1000, 200, 200] (d_e = 7.5um)
% ! Use avg of hist_mean(d) 120 ~ 400 as baseline...    % X No stim as baseline (just avg for all trials)



% Way of sigma_N_E = 20
Parset_i = 1;
sigma_micron = [144, 146, 110; 3.75, 20, 9.25];
kappa = [0.075, 0.125, 0.125];
K_in = repmat([1000 200 200], 2, 1);
alpha_slow = 0.3;
rX_off = 7.5e-3;    % kHz
Parameters_Recurrent(Parset_i, sigma_micron, kappa, K_in, alpha_slow, rX_off);


% XX
% Way of sigma_N_E = 15 (worse)
Parset_i = 2;
sigma_micron = [144, 146, 110; 3.75, 15, 9.25];
kappa = [0.075, 0.125, 0.125];
K_in = repmat([1000 200 200], 2, 1);
alpha_slow = 0.3;
rX_off = 7.5e-3;    % kHz
Parameters_Recurrent(Parset_i, sigma_micron, kappa, K_in, alpha_slow, rX_off);


% kappa0 control
Parset_i = 3;
sigma_micron = [144, 146, 110; 144, 146, 110];
kappa = [0 0 0];
K_in = repmat([1000 200 200], 2, 1);
alpha_slow = 0.3;
rX_off = 7.5e-3;    % kHz
Parameters_Recurrent(Parset_i, sigma_micron, kappa, K_in, alpha_slow, rX_off);


% Use FtrSpc matrix directly
% But 30000 On trials.
if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
dir_local = [dir_0, '/Study/CompNeuro/Projects/Micro-clustering'];
%
Parset_i = 4;
alpha_slow = 0.3;
rX_off = 7.5e-3;    % kHz
FtrSpcWiring_dir = [dir_local, '/Model_SigCorr/SpkSims/Parameters/Parameters_Recurrent_2.mat'];
Parameters_Recurrent_FtrSpc(Parset_i, alpha_slow, rX_off, FtrSpcWiring_dir);

%% Super unstable so give up.

