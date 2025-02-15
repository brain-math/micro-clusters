function [gOSI, Pref_theta, d_BinCenter, SigCorr_histmean, SigCorr_histsem,...
    par2, CI2, par1, CI1] = FR_Analysis(FR, theta_stim, d_BinStep_micron, d_tot, N_sample)

N = size(FR, 1); N_sqrt = round(sqrt(N)); d0 = d_tot / N_sqrt;

gOSI = abs(sum(bsxfun(@times, FR, exp(2 * sqrt(-1) * theta_stim)), 2) ./ sum(FR, 2));
%
Avg_sin = mean(bsxfun(@times, FR, sin(pi + 2 * theta_stim)), 2);
Avg_cos = mean(bsxfun(@times, FR, cos(pi + 2 * theta_stim)), 2);
Pref_theta = reshape((atan2(Avg_sin, Avg_cos) + pi) / 2, [N_sqrt, N_sqrt]);
%
% Avoid empty neurons for spiking sim.
FR_IdxNonEmpty = all(FR == 0, 2);
FR_IdxNonEmpty = find(FR_IdxNonEmpty == 0)';
idx_sample = FR_IdxNonEmpty(randperm(length(FR_IdxNonEmpty)));
if N_sample <= length(idx_sample)
    idx_sample = sort(idx_sample(1: N_sample));
else
    idx_sample = sort(idx_sample);
    fprintf(['Warning: Only ', num2str(length(idx_sample)),...
        ' non empty neurons, less than ', num2str(N_sample), '.\n']);
end
%
idx_sample_j = ceil(idx_sample / N_sqrt); idx_sample_mi = (N_sqrt + 1) - (mod(idx_sample - 1, N_sqrt) + 1);
% x ~ j, or y ~ mi: d0 * [1: N_sqrt] - d0 / 2;
idx_sample_x = d0 * idx_sample_j - d0 / 2; idx_sample_y = d0 * idx_sample_mi - d0 / 2;
distance = sqrt(bsxfun(@minus, idx_sample_x', idx_sample_x) .^ 2 +...
    bsxfun(@minus,  idx_sample_y', idx_sample_y) .^ 2);
distance = triu_new(distance, 0, 1);
%
%%d_BinEdge = [- d_BinStep_micron / 2: d_BinStep_micron: d_tot + d_BinStep_micron / 2, 5 * d_tot];
%%d_BinCenter = [0: d_BinStep_micron: d_tot]';
%d_BinEdge = [0: d_BinStep_micron: d_tot, 5 * d_tot];
%d_BinCenter = [d_BinStep_micron / 2: d_BinStep_micron: d_tot - d_BinStep_micron / 2]';
if d_BinStep_micron == 0    % Then automatic real dist before 50, and 10 after 50. For narrow.
    d_uniq = sort(unique(distance(distance < 50)));    % columnar vector
    d_BinCenter = [d_uniq; [60: 10: d_tot]'];
    dd = diff(d_BinCenter) / 2;
    d_BinEdge = [d_BinCenter(1) - dd(1); d_BinCenter(1: end - 1) + dd; d_BinCenter(end) + dd(end); 5 * d_tot];
else    % Broad, eg, 15
    d_BinEdge = [- d_BinStep_micron / 2: d_BinStep_micron: d_tot + d_BinStep_micron / 2, 5 * d_tot];
    d_BinCenter = [0: d_BinStep_micron: d_tot]';
end
%
FR_sample = FR(idx_sample', :);
SigCorr = corrcoef(FR_sample');
SigCorr = triu_new(SigCorr, 0, 1);
%
[SigCorr_histmean, SigCorr_histsem, N_hist] = histogram_mean_sem(SigCorr, distance, d_BinEdge);
%
idx_start = min(find(SigCorr_histmean > 0));
SigCorr_histmean = SigCorr_histmean(idx_start: end - 1);
SigCorr_histsem = SigCorr_histsem(idx_start: end - 1);
d_BinCenter = d_BinCenter(idx_start: end);
Empty_idx = find(SigCorr_histmean == 0);
if ~isempty(Empty_idx)
for k = 1: length(Empty_idx)
SigCorr_histmean(Empty_idx(k)) =...
    (SigCorr_histmean(Empty_idx(k) - 1) +...
    SigCorr_histmean(Empty_idx(k) + 1)) / 2;
end
end
%
%
idx_valid = find(d_BinCenter <= 500);    %find(SigCorr_histsem ~= 0);
%
Exp2 = @(x, A, lambda, b) A * exp(- (x / lambda).^2) + b;
Err = @(par) (Exp2(d_BinCenter(idx_valid), par(1), par(2), par(3)) -...
    SigCorr_histmean(idx_valid)) ./ SigCorr_histsem(idx_valid);
IniVal = [1, 10, 0]; par_lb = [1e-5, 1e-5, 0]; par_ub = [1000, 1000, 1000];
options = optimoptions('lsqnonlin', 'Display', 'none', 'MaxFunEvals', 1200, 'MaxIter', 1200);
[par2, ~, residual, exitflag, ~, ~, Jacobian] = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
CI2 = nlparci(par2, residual, 'jacobian', Jacobian);
%
Exp1 = @(x, A, lambda, b) A * exp(- (x / lambda)) + b;
Err = @(par) (Exp1(d_BinCenter(idx_valid), par(1), par(2), par(3)) -...
    SigCorr_histmean(idx_valid)) ./ SigCorr_histsem(idx_valid);
[par1, ~, residual, exitflag, ~, ~, Jacobian] = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
CI1 = nlparci(par1, residual, 'jacobian', Jacobian);

