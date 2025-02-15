function [d_BinCenter, NoiseCorr_histmean, NoiseCorr_histsem,...
    par2, CI2, N_valid] = NoiseCorr_Analysis(FR, d_BinStep_micron, d_tot, FR_thr)

N = size(FR, 1); N_sqrt = round(sqrt(N)); d0 = d_tot / N_sqrt;
N_stim = size(FR, 3);


% Bin
idx_j = ceil([1: N] / N_sqrt); idx_mi = (N_sqrt + 1) - (mod([1: N] - 1, N_sqrt) + 1);
% x ~ j, or y ~ mi: d0 * [1: N_sqrt] - d0 / 2;
idx_x = d0 * idx_j - d0 / 2; idx_y = d0 * idx_mi - d0 / 2;
distance = sqrt(bsxfun(@minus, idx_x', idx_x) .^ 2 +...
    bsxfun(@minus,  idx_y', idx_y) .^ 2);
distance = triu_new(distance, 0, 1);
Nd_max = length(distance);
%
if d_BinStep_micron == 0    % Then automatic real dist before 50, and 10 after 50. For narrow.
    d_uniq = sort(unique(distance(distance < 50)));    % columnar vector
    d_BinCenter = [d_uniq; [60: 10: d_tot]'];
    dd = diff(d_BinCenter) / 2;
    d_BinEdge = [d_BinCenter(1) - dd(1); d_BinCenter(1: end - 1) + dd; d_BinCenter(end) + dd(end); 5 * d_tot];
    clear d_uniq dd
else    % Broad, eg, 15
    d_BinEdge = [- d_BinStep_micron / 2: d_BinStep_micron: d_tot + d_BinStep_micron / 2, 5 * d_tot];
    d_BinCenter = [0: d_BinStep_micron: d_tot]';
end
clear idx_j idx_mi idx_x idx_y distance



%FR_thr = 0;
NoiseCorr_histmean_all = NaN(length(d_BinCenter) + 1, N_stim);
NoiseCorr_histsem_all = NaN(length(d_BinCenter) + 1, N_stim);
N_hist = NaN(length(d_BinCenter) + 1, N_stim);
N_valid = NaN(1, N_stim);

for stim_k = 1: N_stim
    FR_k = FR(:, :, stim_k);
    %
    % Threshold of FR
    idx_valid = find(max(FR_k, [], 2) > FR_thr);    %find(mean(FR_k, 2) > FR_thr);
    N_valid(stim_k) = length(idx_valid);
    %
    % Find distance
    idx_j = ceil(idx_valid / N_sqrt); idx_mi = (N_sqrt + 1) - (mod(idx_valid - 1, N_sqrt) + 1);
    % x ~ j, or y ~ mi: d0 * [1: N_sqrt] - d0 / 2;
    idx_x = d0 * idx_j - d0 / 2; idx_y = d0 * idx_mi - d0 / 2;
    distance = sqrt(bsxfun(@minus, idx_x', idx_x) .^ 2 +...
        bsxfun(@minus,  idx_y', idx_y) .^ 2);
    distance = triu_new(distance, 0, 1);
    %
    NoiseCorr = triu_new(corrcoef(FR_k(idx_valid, :)'), 0, 1);
    [NoiseCorr_histmean_all(:, stim_k), NoiseCorr_histsem_all(:, stim_k), N_hist(:, stim_k)] =...
        histogram_mean_sem(NoiseCorr, distance, d_BinEdge);
    clear FR_k idx_valid idx_j idx_mi idx_x idx_y distance NoiseCorr
end
NoiseCorr_histmean_all = NoiseCorr_histmean_all(1: end - 1, :);
NoiseCorr_histsem_all = NoiseCorr_histsem_all(1: end - 1, :);
N_hist = N_hist(1: end - 1, :);
%
N_tot = sum(N_hist, 2);
NoiseCorr_histmean = sum(NoiseCorr_histmean_all .* N_hist, 2) ./ N_tot;
Std = NoiseCorr_histsem_all .* sqrt(N_hist);
Sum_sample_square = sum((N_hist - 1) .* (Std .^ 2) + N_hist .* (NoiseCorr_histmean_all .^ 2), 2);
NoiseCorr_histsem = sqrt((Sum_sample_square - N_tot .* (NoiseCorr_histmean .^ 2)) ./ (N_tot .* (N_tot - 1)));



idx_start = min(find(NoiseCorr_histmean > 0));
NoiseCorr_histmean = NoiseCorr_histmean(idx_start: end);
NoiseCorr_histsem = NoiseCorr_histsem(idx_start: end);
d_BinCenter = d_BinCenter(idx_start: end);
Empty_idx = find(NoiseCorr_histmean == 0);
if ~isempty(Empty_idx)
for k = 1: length(Empty_idx)
    NoiseCorr_histmean(Empty_idx(k)) =...
        (NoiseCorr_histmean(Empty_idx(k) - 1) +...
        NoiseCorr_histmean(Empty_idx(k) + 1)) / 2;
end
end


idx_valid = find(d_BinCenter <= 500);    %find(SigCorr_histsem ~= 0);
%
DualExp2 = @(x, A1, A2, lambda1, lambda2, b)...
    A1 * exp(- (x / lambda1).^2) + A2 * exp(- (x / lambda2).^2) + b;
IniVal = [0.5, 0.05, 10, 200, 0];
par_lb = [0, 0, 0, 0, -10]; par_ub = [100, 100, 10000, 10000, 10];
options = optimoptions('lsqnonlin', 'Display', 'none',...
    'MaxFunEvals', 1200, 'MaxIter', 1200, 'StepTolerance', 1e-8);
Err = @(par) DualExp2(d_BinCenter(idx_valid), par(1), par(2), par(3), par(4), par(5)) -...
    NoiseCorr_histmean(idx_valid);    % Break if ./NoiseCorr_histsem(idx_valid)
[par2, ~, residual, exitflag, ~, ~, Jacobian] = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
CI2 = nlparci(par2, residual, 'jacobian', Jacobian);


