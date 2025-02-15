function sf = FFWD_spiking_generation_LIF(Nf, rX_k, J_matrix, rF_off, T_on, T_off, T)

%rX_k = rX(:, Feature_i);
sf = zeros(2, round(Nf * T * max(rX_k) * 5));
% The spiking record for F will be saved in a 2-row matrix. 
% The # of column indicates the # of all spikings of all neurons.
% The 1st row contains spike times, and the 2nd row contains indices of neurons that spike.
nspks = 0;
T_duration = T_on + T_off;

V_reset = 0; V_thr = 1;
sigma_wn = 15000 * sqrt(0.1e-3 / 1e-3);    % For dt = 1ms.

for t = 1: T
if mod(t, T_duration) < T_off + 0.1    % OFF duration
    rF_t = ones(Nf, 1) * rF_off;    % in kHz
    ifspk = (rand(Nf, 1) < rF_t);
    % Each t interval represents 1ms. P(spike) = 1(ms) x FR(kHz). So use (0, FR) within (0, 1) to sample this prob.
else    % ON duration.
    if mod(t, T_duration) < T_off + 1.5    % 1st msec in ON, initialization
        V = V_reset * ones(Nf, 1);
    end
    x_wn_t = sigma_wn * randn(Nf, 1) * sqrt(1e-3);
    I2F_t = J_matrix * (rX_k * 1e3 + x_wn_t);    % rX_i input is in kHz
    dVdt = -V + 1.1 * I2F_t;
    % !! 1.1 is specific for rX level and sigma_wn = 15000 * sqrt(0.1e-3 / 1e-3),
      % so that f-I curve is approximately y = x !
    V = V + 1e-3 * dVdt;
    %
    ifspk = (V >= V_thr); V(ifspk) = V_reset;
end
% common part
nspks_temp = nnz(ifspk);
if (nspks_temp >= 1)
    sf(1, nspks + (1: nspks_temp)) = t;
    sf(2, nspks + (1: nspks_temp)) = find(ifspk > 0);
    nspks = nspks + nspks_temp;
end
end

[~, j] = sort(sf(1, :)); sf = sf(:, j);    % sort with time.
sf = sf(:, sf(1, :) > 0 & sf(1, :) <= T);    % cut unused tails.

