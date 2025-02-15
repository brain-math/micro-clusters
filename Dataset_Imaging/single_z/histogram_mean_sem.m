function [Value_target_hist2d_mean, Value_target_hist2d_sem,...
    Value_target_N] = histogram_mean_sem(Value_target, Value_coor, Value_coor_BinEdge)
% Normal 2D histogram gives the NUMBER of variables in Value_target fall in each Value_coor_bin.
% But we want the MEAN, STD, .... of that!
% ! Value_target, Value_coor should be in the same N x 1 size.
    % Value_coor_BinEdge: 1D vector.
%% May fail when bins >> datas and empty bins appear !

if (min(Value_coor_BinEdge) > min(Value_coor)) | (max(Value_coor_BinEdge) < max(Value_coor))
    error('Value coor BinEdge is not broaden enough to cover all Value coor !');
    % if so, 0 may occur in Value_coor_BinIdx, leading to error...
end

if (size(Value_target, 2) ~= 1) || (size(Value_coor, 2) ~= 1)
    error('Value_target, Value_coor should be N x 1.');
end

if size(Value_target, 1) ~= size(Value_coor, 1)
    error('The length of Value_target, Value_coor should be the same. Your dataset may include some empty points.');
end

% Get Value_coor_BinIdx
Value_coor_BinIdx = zeros(size(Value_coor));
for i = 1: (length(Value_coor_BinEdge) - 1)
    Value_coor_lb = Value_coor_BinEdge(i); Value_coor_ub = Value_coor_BinEdge(i + 1);
    Idx_within =  find((Value_coor - Value_coor_lb) .* (Value_coor - Value_coor_ub) < 0);
    Value_coor_BinIdx(Idx_within) = i;
    Idx_within_left_edge = find(Value_coor == Value_coor_lb);
    Value_coor_BinIdx(Idx_within_left_edge) = i;
end
Idx_max_edge = find(Value_coor == max(Value_coor_BinEdge));
Value_coor_BinIdx(Idx_max_edge) = length(Value_coor_BinEdge) - 1;

Value_target_BinIdx = [Value_target, Value_coor_BinIdx];
% Do not find idx in Value_coor_BinIdx one by one... 
% But find "BinEdges" of each distinct values in Value_coor_BinIdx, after sorted.
Value_target_BinIdx = sortrows(Value_target_BinIdx, 2);
edges = find(diff(Value_target_BinIdx(:, 2)) ~= 0);
edges = [0; edges; size(Value_target_BinIdx, 1)];

Value_target_hist2d_mean = zeros(length(Value_coor_BinEdge) - 1, 1);
Value_target_hist2d_sem = zeros(length(Value_coor_BinEdge) - 1, 1);
Value_target_N = zeros(length(Value_coor_BinEdge) - 1, 1);
% You can add any types of statistics here.
for k = 1: (length(edges) - 1)
    idx = Value_target_BinIdx(edges(k) + 1, 2); array_k = Value_target_BinIdx(edges(k) + 1: edges(k + 1), 1);
    Value_target_hist2d_mean(idx) = mean(array_k);
    err = array_k - mean(array_k); Value_target_hist2d_sem(idx) = std(err) / sqrt(length(err));
    Value_target_N(idx) = length(array_k);
end

%Value_coor_BinCenter = zeros(size(Value_coor_BinEdge));
%for i = 1: length(Value_coor_BinEdge) - 1
%    Value_coor_BinCenter(i) = (Value_coor_BinEdge(i) + Value_coor_BinEdge(i + 1)) / 2;
%end
%Value_coor_BinCenter = Value_coor_BinCenter(1: end - 1);

% Then we can: errorbar(Value_coor_BinCenter, Value_target_hist2d_mean', Value_target_hist2d_sem');

