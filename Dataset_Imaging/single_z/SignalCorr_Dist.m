function [SigCorr_1d, dist_1d] = SignalCorr_Dist(FR, Location)

[N_neuron, N_stim] = size(FR);

if size(FR, 1) ~= size(Location, 1)
    error('size(FR, 1) ~= size(Location, 1)');
end

% Get the cotical distance of all neuron pairs.
delta_x = bsxfun(@minus, Location(:, 1), Location(:, 1)');
% A matrix of delta_x for all pairs. size; (N_neuron, N_neuron)
delta_y = bsxfun(@minus, Location(:, 2), Location(:, 2)');
dist_2d = sqrt(delta_x .^ 2 + delta_y .^ 2);
dist_1d = triu_new(dist_2d, 0, 1);

% Get the signal correlation of all neuron pairs.
SigCorr_2d = corrcoef(FR');
SigCorr_1d = triu_new(SigCorr_2d, 0, 1);