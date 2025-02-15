function sx = FFWD_spiking_generation(Nx, rX_off, T)

sx = zeros(2, round(Nx * T * rX_off * 5));
% The spiking record for F will be saved in a 2-row matrix. 
% The # of column indicates the # of all spikings of all neurons.
% The 1st row contains spike times, and the 2nd row contains indices of neurons that spike.
nspks = 0;

for t = 1: T
    ifspk = (rand(Nx, 1) < rX_off);
    % Each t interval represents 1ms. P(spike) = 1(ms) x FR(kHz). So use (0, FR) within (0, 1) to sample this prob.
    nspks_temp = nnz(ifspk);
    if (nspks_temp >= 1)
        sx(1, nspks + (1: nspks_temp)) = t;
        sx(2, nspks + (1: nspks_temp)) = find(ifspk > 0);
        nspks = nspks + nspks_temp;
    end
end

[~, j] = sort(sx(1, :)); sx = sx(:, j);    % sort with time.
sx = sx(:, sx(1, :) > 0 & sx(1, :) <= T);    % cut unused tails.

