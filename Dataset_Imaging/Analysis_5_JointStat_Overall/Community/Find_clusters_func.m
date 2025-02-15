function [Pref_orientation_OS_o, Comm_neuron_idx, Size_comm, idx_within_comm, idx_beyond_comm] =...
    Find_clusters_func(ROI_centeroid_OS_o, Pref_orientation_OS_o, Block_N_OS_o,...
    dist_thr, dtheta_thr, if_shuffle_pref_orientation)
% Of course community must be within the same z -- we're looking for horitonzal clusterings, rather than columns.
% Why we save Pref_orientation_OS_o: in case you want to see shuffled community

N_neuron = length(Pref_orientation_OS_o);
%
if if_shuffle_pref_orientation == 1
    Pref_orientation_OS_o = Pref_orientation_OS_o(randperm(N_neuron));
end
%
neuron_idx_edge = [0 cumsum(Block_N_OS_o)];


comm_k = 0; Comm_neuron_idx = cell(N_neuron, 1);
% Each element contains the neuron index of each nearby community.

for z_k = 1: length(Block_N_OS_o)
neuron_k_idx = (neuron_idx_edge(z_k) + 1): neuron_idx_edge(z_k + 1);
%
ROI_centeroid_OS_k = ROI_centeroid_OS_o(neuron_k_idx, :);
Pref_orientation_OS_k = Pref_orientation_OS_o(neuron_k_idx);
delta_i = bsxfun(@minus, ROI_centeroid_OS_k(:, 1), ROI_centeroid_OS_k(:, 1)');
delta_j = bsxfun(@minus, ROI_centeroid_OS_k(:, 2), ROI_centeroid_OS_k(:, 2)');
dist = sqrt(delta_i .^ 2 + delta_j .^ 2);
dtheta_pref = abs(bsxfun(@dtheta, Pref_orientation_OS_k * (pi / 180),...
    Pref_orientation_OS_k' * (pi / 180))) * (180 / pi);
%
for neuron_k = 1: Block_N_OS_o(z_k)
    idx_nearby = find(dist(neuron_k, :) <= dist_thr);
    idx_nearby_similar_theta = idx_nearby(dtheta_pref(neuron_k, idx_nearby) <= dtheta_thr);
    % For "clusters" detected with more than 2 elements, choose subset whose all subsets have d < d_thr.
    n = length(idx_nearby_similar_theta);
    if n >= 2
        % Get all subsets with >= 2 elements.
        subsets = logical(dec2bin(1: 2^n-1) - '0'); subsets = subsets(sum(subsets, 2) >= 2, :);
        for subset_k = 1: size(subsets, 1)
            idx_subset = idx_nearby_similar_theta(subsets(subset_k, :));
            dist_subset = dist(idx_subset, idx_subset);
            if all(dist_subset <= dist_thr)
                comm_k = comm_k + 1;
                Comm_neuron_idx{comm_k} = idx_subset + neuron_idx_edge(z_k);
            end
        end
    end
end
end
Comm_neuron_idx = Comm_neuron_idx(1: comm_k);%(~cellfun(@isempty, Comm_neuron_idx));

% If someone is the subset of some other one (especially, equal), remove it.
% Sort first.
Comm_neuron_idx_tmp = cell(comm_k, 2);
Comm_neuron_idx_tmp(:, 1) = Comm_neuron_idx;
for k = 1: comm_k, Comm_neuron_idx_tmp{k, 2} = Comm_neuron_idx{k}(1); end
Comm_neuron_idx_tmp = sortrows(Comm_neuron_idx_tmp, 2);
Comm_neuron_idx = Comm_neuron_idx_tmp(:, 1);
idx1 = zeros(1, comm_k);
for k = 1: comm_k, idx1(k) = Comm_neuron_idx_tmp{k, 2}; end
idx1_edge = [0, find(diff(idx1) ~= 0), length(idx1)];
clear Comm_neuron_idx_tmp idx1 comm_k
% clear -- within the same idx(1)
for k = 1: length(idx1_edge) - 1
    % (idx1_edge(k) + 1): idx1_edge(k + 1)
    for i = (idx1_edge(k) + 2): idx1_edge(k + 1)    % to be tested
    for j = (idx1_edge(k) + 1): (i - 1)    % valid standards, before him.
        tmp = intersect(Comm_neuron_idx{i}, Comm_neuron_idx{j});
        if isequal(tmp, Comm_neuron_idx{i})
            Comm_neuron_idx{i} = [];
        elseif (~isequal(tmp, Comm_neuron_idx{i})) && (isequal(tmp, Comm_neuron_idx{j}))
            Comm_neuron_idx{j} = [];
        end
    end
    end
end
Comm_neuron_idx = Comm_neuron_idx(~cellfun(@isempty, Comm_neuron_idx));
% Some >= 3 community may have subsets repeated.
Size_comm = zeros(1, length(Comm_neuron_idx));
for k = 1: length(Comm_neuron_idx)
    Size_comm(k) = length(Comm_neuron_idx{k});
end
%
sizes = unique(Size_comm); sizes = sizes(2: end);    % 3, 4, ...
idx_valid = ones(1, length(Comm_neuron_idx));
%
for size_k = sizes(end: -1: 1)
    idx_big = find(Size_comm == size_k);    % eg. 4
    idx_small = find(Size_comm < size_k);    % eg. 2 and 3. Next round: 3 -> 2.
    for i = 1: length(idx_big)    % Valid, big size standards
    for j = 1: length(idx_small)    % to be tested
        tmp = intersect(Comm_neuron_idx{idx_big(i)}, Comm_neuron_idx{idx_small(j)});
        if isequal(tmp, Comm_neuron_idx{idx_small(j)}), idx_valid(idx_small(j)) = 0; end
    end
    end
end
Comm_neuron_idx = Comm_neuron_idx(find(idx_valid == 1));
Size_comm = Size_comm(find(idx_valid == 1));


% Then findout idx within and without communities, and size of each community.
idx_within_comm = []; 
for k = 1: length(Comm_neuron_idx)
    idx_within_comm = [idx_within_comm, Comm_neuron_idx{k}];
end
idx_within_comm = sort(unique(idx_within_comm));
idx_beyond_comm = sort(setdiff(1: N_neuron, idx_within_comm));
% setdiff(A, B): Find the values in A that are not in B.

