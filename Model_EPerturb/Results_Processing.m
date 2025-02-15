function Results_Processing

Parset_i = str2num(getenv('SLURM_ARRAY_TASK_ID'));
%load([pwd, '/Parameters/Parameters_Recurrent_1.mat'], 'param');
%Simul_idx = (Parset_i - 1) * sum(param.N_simul) + [1: sum(param.N_simul)];
if Parset_i == 1
    Simul_idx = 1: 301;
elseif Parset_i == 2
    Simul_idx = 302: 602;
elseif Parset_i == 3
    Simul_idx = 603: 903;
elseif Parset_i == 4
    Simul_idx = 904: 1504;
end
% sbatch --array=1-301%38 main.sh
% sbatch --array=302-602%22 main.sh
% sbatch --array=603-903%22 main.sh
% sbatch --array=904-1504%22 main.sh
%
load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param');

NTw = floor(param.T / param.Tw);
Tw_per_trial = floor((param.T_off + param.T_on) / param.Tw);    % 6
NTw_valid = NTw - param.N_trial_burn * Tw_per_trial;
FR_Tw_off = NaN(param.N, NTw_valid, param.N_simul(1));
%FR_Tw_on = NaN(param.N, NTw_valid, param.N_simul(2));
for simul_i = 1: param.N_simul(1)
    filename = [pwd, '/Results/Results_Simul', num2str(Simul_idx(simul_i)), '.mat'];
    load(filename, 'FR_Tw'); FR_Tw_off(:, :, simul_i) = FR_Tw;
    %
    if simul_i == 1, load(filename, 'spktrain'); spktrain_off_sample = spktrain; end
end
tmp = FR_Tw_off(1: param.Ne, :, :); FR_E_off = mean(tmp(:)); clear tmp FR_Tw_off
%
FR_E_on = NaN(param.Ne, param.N_trial_per_simul, param.N_simul(2));
for simul_i = 1: param.N_simul(2)
    simul_idx_i = Simul_idx(simul_i + param.N_simul(1));
    load([pwd, '/Results/Results_Simul', num2str(simul_idx_i), '.mat'], 'FR_Tw');
    for trial_i = 1: param.N_trial_per_simul
        Tw_idx_on = (trial_i - 1) * Tw_per_trial + param.Tw_pert;    % [6], or [5, 6]
        FR_E_on(:, trial_i, simul_i) = mean(FR_Tw(1: param.Ne, Tw_idx_on), 2);
    end
    %
    if simul_i == 1, load(filename, 'spktrain'); spktrain_on_sample = spktrain; end
end
FR_E_on = reshape(FR_E_on, param.Ne, []);    % (Ne, N_trial_tot(2))

if Parset_i ~= 4
    save([pwd, '/Results_Total_', num2str(Parset_i), '.mat'],...
        'FR_E_off', 'FR_E_on', 'spktrain_off_sample', 'spktrain_on_sample');
elseif Parset_i == 4
    N_trial_tot = size(FR_E_on, 2);
    FR_E_on_1 = FR_E_on(:, 1: N_trial_tot / 2);
    FR_E_on_2 = FR_E_on(:, (N_trial_tot / 2 + 1): end);
    clear FR_E_on
    save([pwd, '/Results_Total_', num2str(Parset_i), '_A.mat'],...
        'FR_E_off', 'FR_E_on_1', 'spktrain_off_sample', 'spktrain_on_sample');
    save([pwd, '/Results_Total_', num2str(Parset_i), '_B.mat'], 'FR_E_on_2');
end



import java.lang.System
java.lang.System.exit(0)

