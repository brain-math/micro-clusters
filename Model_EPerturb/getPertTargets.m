function [tPert, PertIdx] = getPertTargets(T, T_off, T_on, N_spk_pert, dt, Ne, N_simul_tot)


% Row vector, containing spike times:
% [2000 + 0, 2000 + 100, ..., 2000 + 900, 5000 + 0, ..., 5000 + 900, ...] / 0.05;    (0 idx in C)
% % Also, there should be some variability, or responses of each trial would be identical.
% % i.e. add a random number +0 ~ +10.

N_trial_per_simul_tot = T / (T_off + T_on);

tmp1 = (T_off + T_on) * ([1: N_trial_per_simul_tot] - 1) + T_off;
tmp2 = repmat(tmp1, [N_spk_pert 1]);
interval = T_on / N_spk_pert;
tmp3 = bsxfun(@plus, tmp2, [0: interval: T_on - interval]');
spk_timing_const = tmp3(:)';
tPert = spk_timing_const / dt;
clear tmp1 tmp2 tmp3

%timing_error = (randi(11, [1, N_spk_pert * N_trial_per_simul_tot]) - 1);
%% randi(11, size) - 1 generates random numbers among [0, 1, ..., 10].
%t_pert = (spk_timing_const + timing_error) / dt;
%% so it would be int of ti in the timeloop.

tmp1 = randi(Ne, [N_simul_tot, N_trial_per_simul_tot]);
PertIdx = NaN(N_simul_tot, N_trial_per_simul_tot * N_spk_pert);
for k = 1: N_simul_tot
    tmp2 = repmat(tmp1(k, :), [N_spk_pert 1]);
    PertIdx(k, :) = tmp2(:)';
end
clear tmp1 tmp2
