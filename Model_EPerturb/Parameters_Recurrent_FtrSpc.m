function Parameters_Recurrent(Parset_i, alpha_slow, rX_off, mat_dir)

load(mat_dir, 'param', 'Wrf', 'Wrr');%, 'sigma_micron', 'kappa', 'ifDthetaDep', 'sigma_dtheta', 'a_dtheta');
param_f = param; clear param
Wrx = Wrf; clear Wrf


%% Common paramters of L2/3 Exponential - Integrate - Fire dynamics (except connectivity)
Nx_sqrt = 100; Ne_sqrt = Nx_sqrt; Ni_sqrt = Ne_sqrt / 2;
param.Nx = Nx_sqrt ^ 2; param.Ne = Ne_sqrt ^ 2; param.Ni = Ni_sqrt ^ 2; param.N = param.Ne + param.Ni;
param.d_tot = Nx_sqrt * 7.5;
%
param.N_trial_tot = [50, 30000];    % More !
param.N_trial_per_simul = 50;    % ~ 1h
param.N_simul = param.N_trial_tot / param.N_trial_per_simul;
param.N_trial_burn = 1;
param.T_off = 2000;
param.T_on = 1000;
param.T = (param.N_trial_per_simul + param.N_trial_burn) * (param.T_off + param.T_on);
param.dt = 0.05;
%
param.N_spk_pert = 10;    % 10 spikes within T_on 1s, i.e. 10 Hz
param.Tw = 500;
param.Tw_pert = [6];    % [5 6]
[param.tPert, param.PertIdx_Tot] = getPertTargets(param.T, param.T_off, param.T_on,...
    param.N_spk_pert, param.dt, param.Ne, sum(param.N_simul));
%param.PertIdx = ...
%
param.Cm = [1 1];
param.gl = [1/15 1/10];    % g_leak. Actually g_leak_real / C_m. Unit: 1/ms.
param.vl = [-60 -60];    % E_L. Unit: mV.
param.DeltaT = [2 0.5];    % Unit: mV.
param.vT = [-50 -50];    % Unit: mV.
param.vth = [-10 -10];    % Fire threshold. Unit: mV.
param.vre = [-65 -65];    % Unit: mV.
param.tref = [1.5 0.5];    % Unit: ms.
param.vlb = [-100 -100];    % Lower bound of v(t):
% v[j]+=fmax(  ( Isyntot-gl[1]*(v[j]-Vleak[1]) + gl[1]*DeltaT[1]*EXP((v[j]-VT[1])/DeltaT[1]) ) *dt/C[1], Vlb[1]-v[j]  );
%
% Time constants in synaptic temporal kerners \eta_F(t), \eta_E(t), \eta_I(t)
param.taursyn = [1, 10;
                 1, 10;
                 1, 10];
% Synaptic decay time const: \tau_Fd, \tau_Ed, \tau_Id
param.taudsyn = [5, 100;
                 5, 100;
                 8, 100];
% Fast / slow combination, percentage of fast
param.alpha_slow = alpha_slow;
param.Psyn = [1 - param.alpha_slow, param.alpha_slow;
              1 - param.alpha_slow, param.alpha_slow;
              1 - param.alpha_slow, param.alpha_slow];
%
% "maxns": The max possible # of all spikings of all neurons. For the initialization of s1 in C.
param.maxns = param.N * param.T * 0.05;    % overall avg 0.05kHz / 50Hz.
%
% Neuron indice to record synaptic currents: ALL.
param.Irecord = 1: param.N;


% Poisson L4
param.rX_off = rX_off;    % 4 * 1e-3;    % kHz


%% L2/3 Connectivity
param.W_est = param_f.W_est;
param.G_est = param_f.G_est;
param.K_in = param_f.K_in;
param.Jx = param_f.Jx; param.Jr = param_f.Jr;
param.sigma_micron = param_f.sigma_micron;
param.kappa = param_f.kappa;
param.ifDthetaDep = param_f.ifDthetaDep;
param.sigma_dtheta = param_f.sigma_dtheta;
param.a_dtheta = param_f.a_dtheta;

param.Kx = param_f.Kx;
param.KxMax = param_f.KxMax;
param.KxCumsum = param_f.KxCumsum;
param.Kr = param_f.Kr;
param.KeMax = param_f.KeMax;
param.KiMax = param_f.KiMax;


save([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param', 'Wrx', 'Wrr');
fprintf(['Parameters_Recurrent_', num2str(Parset_i), '.mat saved.\n']);


