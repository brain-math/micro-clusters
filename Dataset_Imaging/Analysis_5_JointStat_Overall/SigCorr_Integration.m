function [SigCorr_1d, dist_1d] = SigCorr_Integration(FR, Location, Block_N)

[N_neuron, N_stim] = size(FR);

Block_idx = [0, cumsum(Block_N)];

if size(FR, 1) ~= size(Location, 1)
    error('size(FR, 1) ~= size(Location, 1)');
end
if size(FR, 1) ~= Block_idx(end)
    error('Total number of Block_N does not match FR and Location.');
end

% Get the cotical distance of all neuron pairs.
delta_x = bsxfun(@minus, Location(:, 1), Location(:, 1)');    % A matrix of delta_x for all pairs. size; (N_neuron, N_neuron)
delta_y = bsxfun(@minus, Location(:, 2), Location(:, 2)');
dist_2d = sqrt(delta_x .^ 2 + delta_y .^ 2);
dist_1d = [];
for k = 1: length(Block_N)
    idx = (Block_idx(k) + 1): Block_idx(k + 1);
    dist_block = dist_2d(idx, idx);
    dist_1d = [dist_1d; triu_new(dist_block, 0, 1)];
end

% Get the signal correlation of all neuron pairs.
SigCorr_2d = corrcoef(FR');
SigCorr_1d = [];
for k = 1: length(Block_N)
    idx = (Block_idx(k) + 1): Block_idx(k + 1);
    SigCorr_block = SigCorr_2d(idx, idx);
    SigCorr_1d = [SigCorr_1d; triu_new(SigCorr_block, 0, 1)];
end

% % Get the mean / se values for each point of d_BinCenter
% [SigCorr_mean, SigCorr_se, Num_Bin] = histogram_mean_sem(SigCorr_1d, dist_1d, d_BinEdge);
% % The size of "SigCorr_mean" and "SigCorr_se": (length(d_BinCenter), 1)
