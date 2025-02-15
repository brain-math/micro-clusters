function x_list = hist2list(x_histc)

x_histc = reshape(x_histc, [1 length(x_histc)]);    % row

idx_valid = find(x_histc ~= 0);
N_valid = length(idx_valid);
if N_valid == 0
    error('Empty input vector, r u fooling me?');
end
x_histc_valid = x_histc(idx_valid);

cdf = cumsum(x_histc_valid);
idx_of_idx_valid = ceil(interp1([0, cdf], [1, 1: N_valid], 1: cdf(end), 'linear'));
x_list = idx_valid(idx_of_idx_valid)';    % output is column vec


% x_list = sort(randi(114, [1 1919]));
% x_hist = histcounts(x_list, 0.5: 1: 114.5);
% x_list_2 = hist2list(x_hist);
% tmp = x_list - x_list_2; sum(tmp .^ 2)
