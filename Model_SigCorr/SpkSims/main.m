function main

clc;
tic;

job_idx = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%
% Seed for rand. Should be different for scripts start running at the same time.
rng('shuffle');
seed_offset_simul = randi(floor(intmax/10));
rng(job_idx + seed_offset_simul);


% For any parameter there are 10 stimuli, each ask for 100 trials.
% Each simulation has 50 valid trials, so 2 simulations for each stim. So totally (N_stim) 10 * 2 = 20 simulations.
load([pwd, '/Parameters/Parameters_Recurrent_1.mat'], 'param');
N_SpkCountperSimul = param.N_SpkCountperSimul;    % 50
load([pwd, '/Parameters/Parameters_FFWD.mat'], 'N_stim');    % 10
N_trial = 100; N_sims_AllTrials = round(N_trial / N_SpkCountperSimul);    % 2
N_simulation = N_stim * N_sims_AllTrials;    % 20
%
% There could be multiple parameters, corresponding to different Parameters_Recurrent_x.mat.
% sbatch --array=1-20% main.sh    % corresponding to Parset_i = 1;
Parset_i = ceil(job_idx / N_simulation);
% 1, .........................., 1; 2, .........................., 2; ......
Simulation_i = mod(job_idx - 1, N_simulation) + 1;
% 1, ........................., 40; 1, ........................., 40; ......
Feature_i = ceil(Simulation_i / N_sims_AllTrials);
% 1, ..., 1, 2, ..., 2, ......, 10; 1, ..., 1, 2, ..., 2, ......, 10; ......
Repeat_i = mod(Simulation_i - 1, N_sims_AllTrials) + 1;
% 1, ..., 4, 1, ..., 4, ......., 4; 1, ..., 4, 1, ..., 4, ......., 4; ......

load([pwd, '/Parameters/Parameters_FFWD.mat'],...
    'Nf', 'rX', 'J_matrix', 'rF_off', 'sigma_n', 'tau_n');
load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'],...
    'param', 'Wrf', 'Wrr');
T = param.T; Ne = param.Ne; Ni = param.Ni; N = param.N;
T_on = param.T_on; T_off = param.T_off;

fprintf(['Begin simulation, ', num2str(toc / 60, '%.3f'), ' min.\n']);

% Generate spiking trains of L4.
rX_i = rX(:, Feature_i);
sf = FFWD_spiking_generation(Nf, rX_i, J_matrix, rF_off,...
    T_on, T_off, sigma_n, tau_n, T);
clear rX rX_i J_matrix
% sf: 2 rows. The # of column indicates the # of all spikings of all neurons.
% The 1st row contains spike times, and the 2nd row contains indices of neurons that spike.
fprintf(['FFWD spikes generated in ' , num2str(toc / 60, '%.3f'), ' min.\n']);

% Voltage initialization.
V0min = param.vre(1); V0max = param.vT(1);
% -65; -50. Upper / Lower bound of stochastic sampling.
param.V0 = (V0max - V0min) .* rand(N, 1) + V0min;

% Time window length for spike counts & Isyn.
Tw = gcd(T_off, T_on);
% 100ms. The largest possible value. (Larger the time window, less memory stress)
NTw = floor(T / Tw);    % 260
%
dt_seg = T_on + T_off;
Nseg = ceil(T / dt_seg);    % 52
Nseg_valid = ceil((T - param.Tburn) / dt_seg);    % 50

% SIMULATION !
[s1, Isyn] = EIF1DRFfastslowSyn_flexK(sf, Wrf, Wrr, param, Tw);
clear Wrf Wrr
% s1: 2 x maxns.
% Isyn: (3 * N) x NTw (260).
Redund = find(s1(2, :) == 0);
if ~isempty(Redund)
    End = Redund(1) - 1; s1 = s1(:, 1: End);    % Cut off elements without any recordings.
    fprintf(['L2/3 simulation finished in ' , num2str(toc / 60, '%.3f'), ' min.\n']);
else
    fprintf('Too many spikes, BOOM!\n');
    error('Too many spikes, BOOM!');
end
save([pwd, '/Results/Results_Par', num2str(Parset_i), '_Feature', num2str(Feature_i),...
    '_Repeat', num2str(Repeat_i), '_tmp.mat'], 's1', 'param',...
    'T_on', 'T_off', 'Ne', 'N', 'Tw', 'NTw', 'Nseg', 'Nseg_valid');


% Spike sample for raster plots
Nseg_raster = 5;
s_sample = s1(:, ((s1(1, :) > param.Tburn) & (s1(1, :) <= (param.Tburn + Nseg_raster * (T_off + T_on)))));
s_sample = sortrows(s_sample', 2);
s_sample(:, 1) = s_sample(:, 1) - param.Tburn;
% figure; scatter(s_sample(:, 1), s_sample(:, 2), 2, 'k', 'filled'); axis([0 1000 1 12500]); axis ij;

% Raw spike counts.
s1E = s1(:, s1(2, :) <= Ne); s1I = s1(:, s1(2, :) > Ne); clear s1
F_Spkcount_tmp = spktime2count(sf, 1: Nf, Tw, NTw, 1); clear sf
E_Spkcount_tmp = spktime2count(s1E, 1: Ne, Tw, NTw, 1); clear s1E    % Ne x NTw
I_Spkcount_tmp = spktime2count(s1I, Ne + 1: N, Tw, NTw, 1); clear s1I    % Ni x NTw

% Spike counts & Isyn counts.
Nrecord = length(param.Irecord);    % Should be N.
L4_FR = zeros(Nf, Nseg_valid); 
L23_FR = zeros(N, Nseg_valid);    % will be in Hz.
Isyn_tot = zeros(Nrecord, Nseg_valid, 3);    % N x 50 x 3    % inputs from F, E, I    % will be in V/s.
%
Tw_NumOff = floor(T_off / Tw);
Tw_NumOn = floor(T_on / Tw);
%
for seg_i = 1: Nseg_valid
    seg_k = Nseg - Nseg_valid + seg_i;
    Idx_on = (seg_k - 1) * (Tw_NumOff + Tw_NumOn) + [Tw_NumOff + 1: Tw_NumOff + Tw_NumOn];
    %
    L4_FR(:, seg_i) = sum(F_Spkcount_tmp(:, Idx_on), 2) / (T_on / 1000);    % Hz  
    L23_FR(1: Ne, seg_i) = sum(E_Spkcount_tmp(:, Idx_on), 2) / (T_on / 1000);    % Hz
    L23_FR(Ne + 1: N, seg_i) = sum(I_Spkcount_tmp(:, Idx_on), 2) / (T_on / 1000);
    Isyn_tot(:, seg_i, 1) = mean(Isyn(1: Nrecord, Idx_on), 2);    % V/s
    Isyn_tot(:, seg_i, 2) = mean(Isyn(Nrecord + 1: 2 * Nrecord, Idx_on), 2);
    Isyn_tot(:, seg_i, 3) = mean(Isyn(2 * Nrecord + 1: 3 * Nrecord, Idx_on), 2);
end
Isyn = Isyn_tot; clear Isyn_tot

% Spike counts of each simulation.
save([pwd, '/Results/Results_Par', num2str(Parset_i), '_Feature', num2str(Feature_i),...
    '_Repeat', num2str(Repeat_i), '.mat'], 'L4_FR', 'L23_FR', 'Isyn', 'N_trial', 's_sample');
delete([pwd, '/Results/Results_Par', num2str(Parset_i), '_Feature', num2str(Feature_i),...
    '_Repeat', num2str(Repeat_i), '_tmp.mat']);
fprintf(['Finished and saved! Total Elapsetime = ', num2str(toc / 60, '%.3f'), ' min.\n']);

import java.lang.System
java.lang.System.exit(0)

