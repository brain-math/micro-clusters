function main

clc;
tic;

fprintf('Begin simulation.\n');

job_idx = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%
% Seed for rand. Should be different for scripts start running at the same time.
rng('shuffle');
seed_offset_simul = randi(floor(intmax/10));
rng(job_idx + seed_offset_simul);

load([pwd, '/Parameters/Parameters_Recurrent_1.mat'], 'param');
N_Simul_per_par = param.N_Simul_per_par;    % 20
Parset_i = ceil(job_idx / N_Simul_per_par);
Simul_i = mod(job_idx - 1, N_Simul_per_par) + 1;
Feature_idx = param.Feature_idx;
clear param
% Each simulation has 50 valid trials.
% For each parameter,
  % 1 - 5: 45 degree stimulus in L4 for trial var, totally 250 trials.
  % 6 - 10: 90 degree.
  % 11 - 15: 135 degree.
  % 16 - 20: 180 degree. So totally 1000 trials for evoked.

load([pwd, '/Parameters/Parameters_FFWD.mat'],...
    'Nf', 'rX', 'J_matrix', 'rF_off', 'sigma_n', 'tau_n');
if (Parset_i == 1) | (Parset_i == 2) | (Parset_i == 3) | (Parset_i == 5)
    load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param', 'Wrf', 'Wrr');
elseif (Parset_i == 4) | (Parset_i == 6)
    load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i - 3), '.mat'], 'param', 'Wrf', 'Wrr');
end
Ne = param.Ne; Ni = param.Ni; N = param.N;
T = param.T; T_on = param.T_on; T_off = param.T_off; Tburn = param.Tburn;

% Generate spiking trains of L4.
rX_i = rX(:, Feature_idx(Simul_i));
if (Parset_i == 1) | (Parset_i == 2) | (Parset_i == 3)
    sf = FFWD_spiking_generation(Nf, rX_i, J_matrix, rF_off, T_on, T_off, sigma_n, tau_n, T);
elseif (Parset_i == 4) | (Parset_i == 5) | (Parset_i == 6)
    sf = FFWD_spiking_generation_LIF(Nf, rX_i, J_matrix, rF_off, T_on, T_off, T);
end
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
NTw = floor(T / Tw); NTw_off = floor(T_off / Tw); NTw_on = floor(T_on / Tw);
N_trial = round(T / (T_on + T_off));
N_trial_valid = round((T - Tburn) / (T_on + T_off));

% SIMULATION !
[s1, Isyn] = EIF1DRFfastslowSyn_flexK(sf, Wrf, Wrr, param, Tw);
clear Wrf Wrr
% s1: 2 x maxns.
% Isyn: (3 * N) x NTw
Redund = find(s1(2, :) == 0);
if ~isempty(Redund)
    End = Redund(1) - 1; s1 = s1(:, 1: End);    % Cut off elements without any recordings.
    save([pwd, '/Results/Results_Par', num2str(Parset_i), '_Simul', num2str(Simul_i), '_tmp.mat'],...
        'sf', 's1', 'Isyn', 'param', 'T_on', 'T_off', 'Nf', 'Ne', 'N', 'Tw', 'NTw',...
        'NTw_off', 'NTw_on', 'N_trial', 'N_trial_valid');
    fprintf(['L2/3 simulation finished in ' , num2str(toc / 60, '%.3f'), ' min.\n']);
else
    fprintf('Too many spikes, BOOM!\n');
    save([pwd, '/Results/Results_Par', num2str(Parset_i), '_Simul', num2str(Simul_i), '_tmp.mat'],...
        'sf', 's1', 'Isyn', 'param', 'T_on', 'T_off', 'Nf', 'Ne', 'N', 'Tw', 'NTw',...
        'NTw_off', 'NTw_on', 'N_trial', 'N_trial_valid');
    error('Too many spikes, BOOM!');
end


% Spike sample for raster plots
N_trial_raster = 5;
s_sample = s1(:, ((s1(1, :) > Tburn) & (s1(1, :) <= (Tburn + N_trial_raster * (T_off + T_on)))));
s_sample = sortrows(s_sample', 2);
s_sample(:, 1) = s_sample(:, 1) - Tburn;
% figure; scatter(s_sample(:, 1), s_sample(:, 2), 2, 'k', 'filled'); axis([0 1000 1 12500]); axis ij;

% Raw spike counts.
s1E = s1(:, s1(2, :) <= Ne); s1I = s1(:, s1(2, :) > Ne); clear s1
F_Spkcount_tmp = spktime2count(sf, 1: Nf, Tw, NTw, 1); clear sf 
E_Spkcount_tmp = spktime2count(s1E, 1: Ne, Tw, NTw, 1); clear s1E    % Ne x NTw
I_Spkcount_tmp = spktime2count(s1I, Ne + 1: N, Tw, NTw, 1); clear s1I    % Ni x NTw

% Spike counts & Isyn counts.
Nrecord = length(param.Irecord);    % Should be N.
L4_FR = zeros(Nf, N_trial_valid); 
L23_FR = zeros(N, N_trial_valid);    % will be in Hz.
Isyn_tot = zeros(Nrecord, N_trial_valid, 3);    % inputs from F, E, I    % will be in V/s.
%
for trial_i = 1: N_trial_valid
    trial_k = N_trial - N_trial_valid + trial_i;
    idx_t_on = (trial_k - 1) * (NTw_off + NTw_on) + [NTw_off + 1: NTw_off + NTw_on];
    %
    L4_FR(:, trial_i) = sum(F_Spkcount_tmp(:, idx_t_on), 2) / (T_on / 1000);    % Hz  
    L23_FR(1: Ne, trial_i) = sum(E_Spkcount_tmp(:, idx_t_on), 2) / (T_on / 1000);    % Hz
    L23_FR(Ne + 1: N, trial_i) = sum(I_Spkcount_tmp(:, idx_t_on), 2) / (T_on / 1000);
    Isyn_tot(:, trial_i, 1) = mean(Isyn(1: Nrecord, idx_t_on), 2);    % V/s
    Isyn_tot(:, trial_i, 2) = mean(Isyn(Nrecord + 1: 2 * Nrecord, idx_t_on), 2);
    Isyn_tot(:, trial_i, 3) = mean(Isyn(2 * Nrecord + 1: 3 * Nrecord, idx_t_on), 2);
end
Isyn = Isyn_tot; clear Isyn_tot

% Spike counts of each simulation.
save([pwd, '/Results/Results_Par', num2str(Parset_i), '_Simul', num2str(Simul_i), '.mat'],...
    'L4_FR', 'L23_FR', 'Isyn', 's_sample');
delete([pwd, '/Results/Results_Par', num2str(Parset_i), '_Simul', num2str(Simul_i), '_tmp.mat']);
fprintf(['Finished and saved! Total Elapsetime = ', num2str(toc / 60, '%.3f'), ' min.\n']);

import java.lang.System
java.lang.System.exit(0)

