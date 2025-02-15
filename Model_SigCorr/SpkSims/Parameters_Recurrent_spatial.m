function Parameters_Recurrent_spatial(Parset_i, W_est, K_in, sigma_micron, kappa, alpha_slow)

%% Common paramters of L2/3 Exponential - Integrate - Fire dynamics (except connectivity)
param.N_SpkCountperSimul = 50;
load([pwd, '/Parameters/Parameters_FFWD.mat'], 'Nf_sqrt', 'd_tot');
if alpha_slow == 0
    param.T_off = 100; param.T_on = 400; param.Tburn = 1000;    % ms
else
    param.T_off = 1000; param.T_on = 2000; param.Tburn = 3000;    % ms
end
T = param.N_SpkCountperSimul * (param.T_on + param.T_off) + param.Tburn;
param.T = T;
param.dt = 0.05;    % Simulation timestep, in ms.
%
Ne_sqrt = Nf_sqrt; Ni_sqrt = Ne_sqrt / 2;
Nf = Nf_sqrt ^ 2; Ne = Ne_sqrt ^ 2; Ni = Ni_sqrt ^ 2; N = Ne + Ni;
param.Nx = Nf; param.Ne = Ne; param.Ni = Ni; param.N = N;
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
param.alpha_slow = alpha_slow;
if alpha_slow == 0
    param.taursyn = [1; 1; 1];    % Synaptic rise time const: \tau_Fr, \tau_Er, \tau_Ir
    param.taudsyn = [5; 5; 8];    % Synaptic decay time const: \tau_Fd, \tau_Ed, \tau_Id
    param.Psyn = [1; 1; 1];    % Fast / slow combination.    % Do not consider for now.
else
    param.taursyn = [1, 10;
                     1, 10;
                     1, 10];
    % Synaptic decay time const: \tau_Fd, \tau_Ed, \tau_Id
    param.taudsyn = [5, 100;
                     5, 100;
                     8, 100];
    % Fast / slow combination, percentage of fast
    param.Psyn = [1 - alpha_slow, alpha_slow;
                  1 - alpha_slow, alpha_slow;
                  1 - alpha_slow, alpha_slow];
end
%
% "maxns": The max possible # of all spikings of all neurons. For the initialization of s1 in C.
param.maxns = N * T * 0.05;    % overall avg 0.05kHz / 50Hz.    % 0.02;
%
% Neuron indice to record synaptic currents: ALL.
param.Irecord = 1: N;


%% L2/3 Connectivity
param.W_est = W_est; param.K_in = K_in;
%
param.G_est = [12, 12, 12; 15, 15, 15];
% K_in_ij * J0_ij = Jtot_ij = W_ij / G_ij
J0 = 1e3 * ((W_est ./ param.G_est) ./ K_in);    % (mV)
param.Jx = J0(:, 1); param.Jr = J0(:, 2: 3);


param.sigma_micron = sigma_micron;
param.kappa = kappa;
dir_JN_spatial = [pwd, '/Parameters/JN_', num2str(Parset_i), '.mat'];
getRecGauMixAdjacencyMatrices(d_tot, Nf_sqrt, Ne_sqrt, Ni_sqrt,...
    K_in, sigma_micron, kappa, dir_JN_spatial);
[Wrf, Wrr, param.Kx, param.KxMax, param.KxCumsum,...
    param.Kr, param.KeMax, param.KiMax] =...
    ReformatAdjacencyMatrices(dir_JN_spatial);


save([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param', 'Wrf', 'Wrr');
fprintf(['Parameters_Recurrent_', num2str(Parset_i), '.mat saved.\n']);


