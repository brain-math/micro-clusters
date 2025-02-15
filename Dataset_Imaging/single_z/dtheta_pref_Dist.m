function [dtheta_pref_1d, dist_1d] = dtheta_pref_Dist(Pref_orientation, Location)
% Pref_orientation must be in deg.

[N_neuron, N_stim] = size(Pref_orientation);

if size(Pref_orientation, 1) ~= size(Location, 1)
    error('size(Pref_orientation, 1) ~= size(Location, 1)');
end

% Get the cotical distance of all neuron pairs.
delta_x = bsxfun(@minus, Location(:, 1), Location(:, 1)');
% A matrix of delta_x for all pairs. size; (N_neuron, N_neuron)
delta_y = bsxfun(@minus, Location(:, 2), Location(:, 2)');
dist_2d = sqrt(delta_x .^ 2 + delta_y .^ 2);
dist_1d = triu_new(dist_2d, 0, 1);

% Get those Δθ_pref
Pref_orientation = Pref_orientation * (pi / 180);
dtheta_pref_2d = abs(bsxfun(@dtheta, Pref_orientation, Pref_orientation')) * (180 / pi);
dtheta_pref_1d = triu_new(dtheta_pref_2d, 0, 1);