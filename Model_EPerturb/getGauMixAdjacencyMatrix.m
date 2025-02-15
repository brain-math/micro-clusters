function JN = getGauMixAdjacencyMatrix(N_sqrt_post, N_sqrt_pre, d_pre, kappa, sigma_n, sigma_b, K_in)
N_post = N_sqrt_post ^ 2;
N_pre = N_sqrt_pre ^ 2;

circrandn = @(mu, sigma, min, max, n) mod(round(sigma * randn(n, 1) + mu) - min, max - min + 1) + min;
alpha_n = sigma_n / d_pre; alpha_b = sigma_b / d_pre;
K_n = round(K_in * kappa); K_b = K_in - K_n;

JN = NaN(N_post, N_pre);

for k_post = 1: N_post
    i_post_on_pre = (mod(k_post - 1, N_sqrt_post) + 1) * (N_sqrt_pre / N_sqrt_post);
    j_post_on_pre = ceil(k_post / N_sqrt_post) * (N_sqrt_pre / N_sqrt_post);
    i_pre_n = circrandn(i_post_on_pre, alpha_n, 1, N_sqrt_pre, K_n); 
    j_pre_n = circrandn(j_post_on_pre, alpha_n, 1, N_sqrt_pre, K_n);
    idx_pre_n = (j_pre_n - 1) * N_sqrt_pre + i_pre_n;
    i_pre_b = circrandn(i_post_on_pre, alpha_b, 1, N_sqrt_pre, K_b);
    j_pre_b = circrandn(j_post_on_pre, alpha_b, 1, N_sqrt_pre, K_b);
    idx_pre_b = (j_pre_b - 1) * N_sqrt_pre + i_pre_b;
    JN(k_post, :) = histcounts([idx_pre_n; idx_pre_b], [0.5: 1: N_pre + 0.5]);
end

