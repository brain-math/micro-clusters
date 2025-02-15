function sf = FFWD_spiking_generation(Nf, rX_k, J_matrix, rF_off,...
    T_on, T_off, sigma_n, tau_n, T)

%rX_k = rX(:, Feature_i);
sf = zeros(2, round(Nf * T * max(rX_k) * 5));
% The spiking record for F will be saved in a 2-row matrix. 
% The # of column indicates the # of all spikings of all neurons.
% The 1st row contains spike times, and the 2nd row contains indices of neurons that spike.
nspks = 0;
T_duration = T_on + T_off;

for t = 1: T
    if mod(t, T_duration) < T_off + 0.1    % OFF duration
        rF_t = ones(Nf, 1) * rF_off;
    elseif mod(t, T_duration) < T_off + 1.5    % 1st msec in ON
        noise = sigma_n * randn(Nf, 1) / sqrt(2 * tau_n);  % initial condition.
        rF_t = J_matrix * (rX_k + noise);
    else    % ON duration.
        dnoise = 1 * (1 / tau_n) * (- noise + sigma_n * randn(Nf, 1));
        noise = noise + dnoise;
        rF_t = J_matrix * (rX_k + noise);
    end
    %
    ifspk = (rand(Nf, 1) < rF_t);
    % Each t interval represents 1ms. P(spike) = 1(ms) x FR(kHz). So use (0, FR) within (0, 1) to sample this prob.
    nspks_temp = nnz(ifspk);
    if (nspks_temp >= 1)
        sf(1, nspks + (1: nspks_temp)) = t;
        sf(2, nspks + (1: nspks_temp)) = find(ifspk > 0);
        nspks = nspks + nspks_temp;
    end
end

[~, j] = sort(sf(1, :)); sf = sf(:, j);    % sort with time.
sf = sf(:, sf(1, :) > 0 & sf(1, :) <= T);    % cut unused tails.

