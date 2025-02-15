function getConnGauF(Nf_sqrt, K, sigma_micron, d_tot, pdf_Thr, savename)

Nf = Nf_sqrt ^ 2; J0 = 1 / K; 
alpha = sigma_micron / d_tot;    % normalized.
%
d = 1 / Nf_sqrt; CoorX = d/2: d: (1 - d/2); CoorY = CoorX(end: -1: 1);
%
rndset = rand(K, Nf);
%
J_matrix = sparse(Nf, Nf);

for j = 1: Nf    % Per post F neurons
    j_post = ceil(j / Nf_sqrt); x_post = CoorX(j_post);
    i_post = mod(j - 1, Nf_sqrt) + 1; y_post = CoorY(i_post);
    delta_x = min(abs([CoorX - x_post; CoorX - (x_post - 1); CoorX - (x_post + 1)]), [], 1);
    delta_y = min(abs([CoorY - y_post; CoorY - (y_post - 1); CoorY - (y_post + 1)]), [], 1)';
    dist = sqrt(bsxfun(@plus, delta_x .^ 2, delta_y .^ 2));
    pdf_raw = (1 / (2 * pi * alpha^2)) * exp(- dist .^2 / (2 * alpha^2));
    pdf_raw = pdf_raw / max(pdf_raw(:));    % Normalize max = 1 to meet pdf_Thr.
    pre_idx = pdf_rand(pdf_raw, pdf_Thr, rndset(:, j));
    J_matrix(j, :) = histcounts(pre_idx, 0.5: 1: Nf + 0.5);
end
J_matrix = J_matrix * J0;

save(savename, 'J_matrix');
% sum(J_matrix, 2);
