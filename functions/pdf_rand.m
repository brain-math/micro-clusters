function idx = pdf_rand(pdf_raw, Thr, rndset)
% rndset = rand(K, 1);

N_col = size(pdf_raw, 1);
K = length(rndset);

[valid_i, valid_j] = find(pdf_raw >= Thr);
N_valid = length(valid_i);
valid_idx_col = (valid_j - 1) * N_col + valid_i;
%
valid_row = sortrows([valid_i, valid_j], 1);
valid_i = valid_row(:, 1); valid_j = valid_row(:, 2); 
valid_idx_row = (valid_j - 1) * N_col + valid_i;

pdf_1d = pdf_raw(valid_idx_col);
cdf = cumsum(pdf_1d); cdf = cdf / cdf(end);
idx_1_raw = round(interp1(cdf, 1: N_valid, rndset(1: round(K/2)), 'linear', 0));
idx_1_raw(idx_1_raw <= 1) = 1; idx_1_raw(idx_1_raw >= N_valid) = N_valid;
idx_1 = valid_idx_col(idx_1_raw);

pdf_1d = pdf_raw(valid_idx_row);
cdf = cumsum(pdf_1d); cdf = cdf / cdf(end);
idx_2_raw = round(interp1(cdf, 1: N_valid, rndset((round(K/2) + 1): end), 'linear', 0));
idx_2_raw(idx_2_raw <= 1) = 1; idx_2_raw(idx_2_raw >= N_valid) = N_valid;
idx_2 = valid_idx_row(idx_2_raw);

idx = sort([idx_1; idx_2]);
