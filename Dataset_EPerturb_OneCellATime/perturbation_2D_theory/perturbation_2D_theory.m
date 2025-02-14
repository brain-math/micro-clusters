if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
addpath(genpath([dir_0, '/Study/CompNeuro/Projects/Functions_simul/']));

% 1) Integration rather than summation, when allocate a constant summation of W into neuron population!
% 2) Ways of efficient and correct convolution!

sigma_b_e = 146; sigma_b_i = 110;
%
% G = diag([12 15]);    % V^(-1)
% J_indiv = [90, -350; 70, -250] / 1000;    % V
% P = [0.01, 0.04; 0.03, 0.04];
% Ne = 10000; Ni = 2500; N = diag([Ne Ni]); 
% W0 = G * (J_indiv .* P) * N / sqrt(Ne + Ni)
%     0.9660   -3.7566
%     2.8174   -3.3541
%Wee0 = 1.08; Wie0 = 2.82; Wei0 = -3.76;
Wee0 = 1.08; Wie0 = 1.25; Wei0 = -3;
kappa_e = 0.05; sigma_n_e = 10; 
kappa_i = 0.05; sigma_n_i = 10;


dxy = 5; xymax = 500; r1 = -xymax: dxy: xymax; idx = find(r1 == 0);
[x, y] = meshgrid(r1); y = y(end: -1: 1, :); r = sqrt(x .^ 2 + y .^ 2);
Gau2 = @(r, sigma) (1 / (2 * pi * sigma^2)) * exp(- r.^ 2 / (2 * sigma^2));    % 2d model
Gau2Mix = @(r, kappa, sigma_n, sigma_b) kappa * Gau2(r, sigma_n) + (1 - kappa) * Gau2(r, sigma_b);
% if you're trying to allocate a constant summation of W into neuron population,
  % Note that all "spatial functions", or xxx ~ distance, represents p.d.f. rather than real probability,
  % so you need integration rather than summation, i.e. * dx / * dx * dy / * dr * (2 * pi * r),
  % or results will be dx/dy dependent -- which must be wrong!
WeeR = Wee0 * Gau2Mix(r, kappa_e, sigma_n_e, sigma_b_e);
WieR = Wie0 * Gau2Mix(r, kappa_e, sigma_n_e, sigma_b_e);
WeiR = Wei0 * Gau2Mix(r, kappa_i, sigma_n_i, sigma_b_i);
%
dr_EE = WeeR;
dr_EEE = conv2_periodic_by_fft2(WeeR, WeeR) * (dxy ^ 2);
% conv2_periodic_by_fft2(WeeR, WeeR) is exactly the same as
  % conv2(WeeR, WeeR, 'same'), but not silly ifft2(fft2(WeeR) .* fft2(WeeR));
dr_EIE = conv2_periodic_by_fft2(WeiR, WieR) * (dxy ^ 2);
drE = dr_EE + dr_EEE + dr_EIE;


% Theory (approx.)
r1p = 0: dxy: xymax;
A = - Wei0 * Wie0 / Wee0;
dr_ee_narrow = Wee0 * kappa_e * Gau2(r1p, sigma_n_e);
dr_ee_broad = Wee0 * Gau2(r1p, sigma_b_e);
dr_ee_Okappa = - Wee0 * kappa_e * Gau2(r1p, sigma_b_e);
%
dr_eee = (Wee0 ^ 2) * Gau2(r1p, sqrt(2) * sigma_b_e);
dr_eee_Okappa = 2 * kappa_e * (Wee0 ^ 2) * (Gau2(r1p, sigma_b_e) - Gau2(r1p, sqrt(2) * sigma_b_e));
dr_eee_Okappa2 = (kappa_e ^ 2) * (Wee0 ^ 2) * Gau2(r1p, sqrt(2) * sigma_n_e);
%
dr_eie = - Wee0 * A * Gau2(r1p, sqrt(sigma_b_e^2 + sigma_b_i^2));
dr_eie_Okappa = - Wee0 * A * kappa_e * (Gau2(r1p, sigma_b_i) - Gau2(r1p, sqrt(sigma_b_e^2 + sigma_b_i^2)))...
    - Wee0 * A * kappa_i * (Gau2(r1p, sigma_b_e) - Gau2(r1p, sqrt(sigma_b_e^2 + sigma_b_i^2)));
dr_eie_Okappa2 = - Wee0 * A * kappa_e * kappa_i * Gau2(r1p, sqrt(sigma_n_e^2 + sigma_n_i^2));
%
dr_theory_O1 = dr_ee_narrow + dr_ee_broad + dr_eee + dr_eie;
dr_theory_Okappa = dr_theory_O1 + dr_ee_Okappa + dr_eee_Okappa + dr_eie_Okappa;
dr_theory_Okappa2 = dr_theory_O1 + dr_eee_Okappa2 + dr_eie_Okappa2;
dr_theory_full = dr_theory_Okappa + dr_eee_Okappa2 + dr_eie_Okappa2;



figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.65, 0.65]);
%
subplot(1, 2, 1); hold on; l = zeros(1, 8);
l(1) = plot(r1(idx: end), dr_EE(idx, idx: end), 'r');
l(2) = plot(r1(idx: end), dr_EEE(idx, idx: end), 'm');
l(3) = plot(r1(idx: end), dr_EIE(idx, idx: end), 'b');
l(4) = plot(r1(idx: end), drE(idx, idx: end), 'k');
l(5) = plot(r1p, dr_theory_O1, 'color', [0 1 0], 'LineStyle', '--');
l(6) = plot(r1p, dr_theory_Okappa, 'color', [0 0.5 0], 'LineStyle', '--');
l(7) = plot(r1p, dr_theory_Okappa2, 'color', [1 0.5 0.5], 'LineStyle', '--');
l(8) = plot(r1p, dr_theory_full, 'k--');
%plot([0 xymax], [0 0], 'k--');
legend(l, {'E - E', 'E - E - E', 'E - I - E', 'Total', 'Theory (O(1) except narrow)',...
    'Theory (+ O(\kappa))', 'Theory (+ O(\kappa^2))', 'Theory (full)'});
xlim([0 xymax]); axis square; grid on;
set(gca, 'XTick', [0: 10: 50, 100: 100: xymax],...
    'XTickLabel', {'0', '', '', '', '', '50', '100', '200', '300', '400', '500'});
YLim1 = get(gca, 'ylim'); YTick1 = get(gca, 'YTick');
xlabel('Distance (\mum)'); ylabel('\DeltaFR');
title('"Simulation" v.s. theory', 'FontWeight', 'normal');
%
subplot(1, 2, 2); hold on; l = zeros(1, 9);
l(1) = plot(r1p, dr_ee_narrow, 'r');
l(2) = plot(r1p, dr_ee_broad, 'color', [0.5 0 0]);
l(3) = plot(r1p, dr_ee_Okappa + 5e-7, 'r--');
l(4) = plot(r1p, dr_eee, 'm');
l(5) = plot(r1p, dr_eee_Okappa, 'm--');
l(6) = plot(r1p, dr_eee_Okappa2, 'color', 'm', 'linestyle', '--');
l(7) = plot(r1p, dr_eie, 'b');
l(8) = plot(r1p, dr_eie_Okappa, 'b--');
l(9) = plot(r1p, dr_eie_Okappa2, 'b--');
%plot([0 xymax], [0 0], 'k--');
legend(l, {'E - E, Narrow: W_{EE}^0 \kappa_E G(\sigma_N^E)',...
     'E - E, O(1) Broad: W_{EE}^0 G(\sigma_B^E)',...
     'E - E, O(\kappa)', ...
     'E - E - E, O(1): W_{EE}^{0 2} G(\surd2 \sigma_B^E)',...
     'E - E - E, O(\kappa)', 'E - E - E, O(\kappa^2)', ...
     'E - I - E, O(1): - W_{EE}^0 A G(\surd(\sigma_B^{E 2} + \sigma_B^{I 2})))',...
     'E - I - E, O(A\kappa)', 'E - I - E, O(A\kappa^2)'});
xlim([0 xymax]); ylim(YLim1); axis square; grid on;
set(gca, 'XTick', [0: 10: 50, 100: 100: xymax],...
    'XTickLabel', {'0', '', '', '', '', '50', '100', '200', '300', '400', '500'}, 'YTick', YTick1);
xlabel('Distance (\mum)'); ylabel('\DeltaFR');
suptitle({['W_{EE}^0 = ', num2str(Wee0), ', W_{IE}^0 = ', num2str(Wie0),...
    ', W_{EI}^0 = ', num2str(Wei0), ',  A = - W_{IE}^0 * W_{EI}^0 / W_{EE}^0 = ', num2str(A, '%.2f')],...
    ['\kappa_E = ', num2str(kappa_e), ', \kappa_I = ', num2str(kappa_i), '; \sigma_N^E = ',...
    num2str(sigma_n_e), ', \sigma_N^I = ', num2str(sigma_n_i), ' (\mum); \sigma_B^E = ',...
    num2str(sigma_b_e), ', \sigma_B^I = ', num2str(sigma_b_i), ' (\mum)']}, 1, 0.95);
%
pause(2); print(gcf, '-dpng', '1.png');
close;

