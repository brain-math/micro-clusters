function [FR_n_L23E, FR_n_L23I, Cov_r_L23E, Cov_r_L23I] =...
    LR_continuum_FT_func(d_micron, fitpar_exp2_input, W, sigma_micron, kappa)

% FT_Hankel = @(r, n, fr) sum(...
%     bsxfun(@times,...
%         besselj(0, 2 * pi * bsxfun(@times, n', r)),...
%         fr .* (2 * pi * r) * (r(2) - r(1))...
%         ),...
%     2)';
% % r, n, and fr all row vectors.

iFT_Hankel = @(r, n, fn) sum(...
    bsxfun(@times,...
        besselj(0, 2 * pi * bsxfun(@times, r', n)),...
        fn .* (2 * pi * n) * (n(2) - n(1))...
        ),...
    2)';
% r, n, and fn all row vectors.

% Such Hankel transformation could be numerically unstable!
% For \int{ J0(2 * pi * n * r) * f(r) * (2 * pi * r) * dr},
  % If there is a small basal value e0 in f(r),
  % You'll have a term e0 * \int{ J0(2 * pi * n * r) * (2 * pi * r) * dr},
  % and this is an amplifying oscillatory term!
% So, fr or fn send to Hankel transformation MUST STRICTLY CONVERGE TO 0!
% As well, make dx as small as possible!
% Last, double check any Hankel transformation steps !!!  (Then why should I use this bullshit?)
  % Or you'll see weird errors......

% We'll use analytical results for L4 part, so L4 corr function could strictly converge to 0 with b = 0.
% Corr_r = A * exp(- (r / lambda) .^ 2);
% -- FT --> Corr_n = (A * pi * lambda ^ 2) * exp(- (pi * lambda * n) .^ 2);
% -- sqrt --> FR_n = sqrt(A * pi * lambda ^ 2) * exp(- (pi * lambda * n) .^ 2 / 2);
% -- iFT --> FR_r = 2 * sqrt(A / (pi * lambda ^ 2)) * exp(- (r / (lambda / sqrt(2))) .^ 2);
% See validation in Exp2FT.nb
%
n_cycmicron = linspace(0, 0.2, 5001);    
FR_n_input = sqrt(fitpar_exp2_input(1) * pi * fitpar_exp2_input(2) ^ 2) *...
    exp(- (pi * fitpar_exp2_input(2) * n_cycmicron) .^ 2 / 2);

% In case you want to check Hankel transformation for L4 part:
%
% Exp2 = @(x, A, lambda, b) A * exp(- (x / lambda).^2) + b;
% Corr_r_L4_0 = Exp2(d_micron, fitpar_exp2_input(1), fitpar_exp2_input(2), 0);
% % You must remove basal value, or collapse !!!
%
% Corr_n_L4 = FT_Hankel(d_micron, n_cycmicron, Corr_r_L4_0);
% % Corr_n_L4_ana = (fitpar_exp2_input(1) * pi * fitpar_exp2_input(2) ^ 2) *...
% %     exp(- (pi * fitpar_exp2_input(2) * n_cycmicron) .^ 2);
% % figure; hold on; plot(n_cycmicron, Corr_n_L4); plot(n_cycmicron, Corr_n_L4_ana, '--');
%
% n_thr_idx = min(find(Corr_n_L4 < Corr_n_L4(1) * 1e-4)); Corr_n_L4(n_thr_idx: end) = 0;
% FR_n_input = sqrt(Corr_n_L4);    % Everything is real, so no worry
% % FR_r_L4 = iFT_Hankel(d_micron, n_cycmicron, FR_n_input);
% % FR_r_L4_ana = 2 * sqrt(fitpar_exp2_input(1) / (pi * fitpar_exp2_input(2)^2)) *...
% %     exp(- 2 * (d_micron / fitpar_exp2_input(2)).^2);
% % figure; hold on; plot(d_micron, FR_r_L4); plot(d_micron, FR_r_L4_ana, '--');
% % n_cycmicron must dense enough or FR_r_L4 will collapse...

Gau2n = @(n, sigma) exp(- (2 * pi * sigma * n) .^ 2 / 2);
Gau2Mixn = @(n, kappa, sigma_n, sigma_b) kappa * Gau2n(n, sigma_n) + (1 - kappa) * Gau2n(n, sigma_b);
%
Wef_n = W(1, 1) * Gau2Mixn(n_cycmicron, kappa(1), sigma_micron(2, 1), sigma_micron(1, 1));
Wee_n = W(1, 2) * Gau2Mixn(n_cycmicron, kappa(2), sigma_micron(2, 2), sigma_micron(1, 2));
Wei_n = W(1, 3) * Gau2Mixn(n_cycmicron, kappa(3), sigma_micron(2, 3), sigma_micron(1, 3));
Wif_n = W(2, 1) * Gau2Mixn(n_cycmicron, kappa(1), sigma_micron(2, 1), sigma_micron(1, 1));
Wie_n = W(2, 2) * Gau2Mixn(n_cycmicron, kappa(2), sigma_micron(2, 2), sigma_micron(1, 2));
Wii_n = W(2, 3) * Gau2Mixn(n_cycmicron, kappa(3), sigma_micron(2, 3), sigma_micron(1, 3));
%
Det = (1 - Wee_n) .* (1 - Wii_n) - Wei_n .* Wie_n;
FR_n_L23E = (1 ./ Det) .* ((1 - Wii_n) .* Wef_n + Wei_n .* Wif_n) .* FR_n_input;
FR_n_L23I = (1 ./ Det) .* (Wie_n .* Wef_n + (1 - Wee_n) .* Wif_n) .* FR_n_input;

Cov_r_L23E = iFT_Hankel(d_micron, n_cycmicron, FR_n_L23E .^ 2);
Cov_r_L23I = iFT_Hankel(d_micron, n_cycmicron, FR_n_L23I .^ 2);

