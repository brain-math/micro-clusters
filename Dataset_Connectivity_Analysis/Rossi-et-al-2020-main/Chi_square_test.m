function p_value = Chi_square_test(y)
% https://online.stat.psu.edu/stat500/lesson/8/8.1

% y = [106, 124, 1;...
%     67, 85, 1;...
%     37, 72, 3];

idx_nonempty_row = ~all(y == 0, 2);
y = y(idx_nonempty_row, :);
% No empty lines allowed in y_null of chi-square test.
DoF = (size(y, 1) - 1) * (size(y, 2) - 1);
y_null = bsxfun(@times, sum(y, 1), sum(y, 2)) / sum(y(:));    % H0: independent
tmp = ((y - y_null) .^ 2) ./ y_null;
chi_square = sum(tmp(:));
p_value = 1 - chi2cdf(chi_square, DoF);    % p_value of H0
