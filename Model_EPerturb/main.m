function main

clc;
tic;

Simul_i = str2num(getenv('SLURM_ARRAY_TASK_ID'));    % Period: (15000 + 50) / 50 = 301, or (30000 + 50) / 50 = 601.
%
% load([pwd, '/Parameters/Parameters_Recurrent_1.mat'], 'param');
% Parset_i = ceil(Simul_i / sum(param.N_simul));
% If_spk_pert = (mod(Simul_i - 1, sum(param.N_simul)) + 1 > param.N_simul(1));
% % first 1 off, then on.
%
% sbatch --array=1-301%38 main.sh
% sbatch --array=302-602%22 main.sh
% sbatch --array=603-903%22 main.sh
% sbatch --array=904-1504%22 main.sh
if (Simul_i >= 1) & (Simul_i <= 301)
    Parset_i = 1;
    load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param');
    if Simul_i > param.N_simul(1), If_spk_pert = 1; else, If_spk_pert = 0; end
    clear param
    Simul_i_mod_par = Simul_i;
elseif (Simul_i >= 302) & (Simul_i <= 602)
    Parset_i = 2;
    load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param');
    if Simul_i > 301 + param.N_simul(1), If_spk_pert = 1; else, If_spk_pert = 0; end
    clear param
    Simul_i_mod_par = Simul_i - 301;
elseif (Simul_i >= 603) & (Simul_i <= 903)
    Parset_i = 3;
    load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param');
    if Simul_i > 602 + param.N_simul(1), If_spk_pert = 1; else, If_spk_pert = 0; end
    clear param
    Simul_i_mod_par = Simul_i - 602;
elseif (Simul_i >= 904) & (Simul_i <= 1504)
    Parset_i = 4;
    load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param');
    if Simul_i > 903 + param.N_simul(1), If_spk_pert = 1; else, If_spk_pert = 0; end
    clear param
    Simul_i_mod_par = Simul_i - 903;
end


rng('shuffle');
seed_offset_simul = randi(floor(intmax/10));
rng(Simul_i + seed_offset_simul);


load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param', 'Wrx', 'Wrr');
%
sx = FFWD_spiking_generation(param.Nx, param.rX_off, param.T);
% sx: 2 rows. The # of column indicates the # of all spikings of all neurons.
% The 1st row contains spike times, and the 2nd row contains indices of neurons that spike.
V0min = param.vre(1); V0max = param.vT(1);
% -65; -50. Upper / Lower bound of stochastic sampling.
param.V0 = (V0max - V0min) .* rand(param.N, 1) + V0min;


% SIMULATION !
fprintf(['L2/3 simulation begin at ' , num2str(toc / 60, '%.3f'), ' min.\n']);
if If_spk_pert == 0
    [spktrain, ~] = EIF1DRFfastslowSyn_flexK(sx, Wrx, Wrr, param, param.Tw);
    % The last input term param.Tw is only used for Isyn. Useless here.
elseif If_spk_pert == 1
    %Simul_i_mod_par = mod(Simul_i - 1, sum(param.N_simul)) + 1;
    param.PertIdx = param.PertIdx_Tot(Simul_i_mod_par, :);
    [spktrain, ~] = EIF1DRFfastslowSyn_flexK_PertE(sx, Wrx, Wrr, param, param.Tw);
end
clear sx Wrx Wrr
Redund = find(spktrain(2, :) == 0);
if ~isempty(Redund)
    End = Redund(1) - 1; spktrain = spktrain(:, 1: End);    % Cut off elements without any recordings.
    fprintf(['L2/3 simulation finished in ' , num2str(toc / 60, '%.3f'), ' min.\n']);
else
    fprintf('Too many spikes, BOOM!\n');
    error('Too many spikes, BOOM!');
end


% Spike counts
Tw = param.Tw;    % 500
NTw = floor(param.T / Tw);
FR_Tw = spktime2count(spktrain, 1: param.N, Tw, NTw, 1) / (Tw / 1000);    % N x NTw, in Hz
NTw_invalid = param.N_trial_burn * (param.T_off + param.T_on) / Tw;
FR_Tw = FR_Tw(:, NTw_invalid + 1: end);
save([pwd, '/Results/Results_Simul', num2str(Simul_i), '.mat'], 'spktrain', 'FR_Tw');
fprintf(['Finished and saved! Total Elapsetime = ', num2str(toc / 60, '%.3f'), ' min.\n']);


import java.lang.System
java.lang.System.exit(0)

