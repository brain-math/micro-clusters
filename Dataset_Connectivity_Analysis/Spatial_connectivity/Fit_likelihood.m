z_boundary = -340;
Pre_Process(z_boundary);
load(['Presynaptic_coordinates_processed_', num2str(-z_boundary), '.mat'],...
    'Pre_L4E', 'Pre_L23E', 'Pre_L23I');

DoubleGau_2 = @(r, kappa, sigma_n, sigma_b)...
    (kappa / (2 * pi * sigma_n^2)) * exp(- r.^2 / (2 * sigma_n^2)) +...
    ((1 - kappa) / (2 * pi * sigma_b^2)) * exp(- r.^2 / (2 * sigma_b^2));
DoubleGau_1 = @(r, kappa, sigma_n, sigma_b)...
    DoubleGau_2(r, kappa, sigma_n, sigma_b) .* (2 * pi * r);
%
opt = statset('Display', 'off', 'MaxIter', 1800, 'MaxFunEvals', 1800, 'TolX', 1e-8, 'TolFun', 1e-8, 'TolBnd', 1e-8);
par_ic = [0.005, 10, 150; 0.01, 30, 150; 0.01, 20, 100];
%
namelist = {'L4E', 'L23E', 'L23I'}; N_type = length(namelist);
par = zeros(N_type, 3); CI = zeros(N_type, 3, 2); %single_sigma = zeros(1, N_type);
for k = 1: N_type
    eval(['r = sqrt(Pre_', namelist{k}, '(:, 1) .^ 2 + Pre_', namelist{k}, '(:, 2) .^ 2);']);
    [par_tmp, CI_tmp] = mle(r, 'pdf', @(r, kappa, sigma_n, sigma_b)DoubleGau_1(r, kappa, sigma_n, sigma_b),...
        'start', par_ic(k, :), 'LowerBound', [0.0001 0 0], 'UpperBound', [0.2, 50, 500], 'options', opt);
    par(k, :) = par_tmp; CI(k, :, 1) = CI_tmp(1, :); CI(k, :, 2) = CI_tmp(2, :);
    %single_sigma(k) = sqrt(mean(r .^ 2) / 2);
end
clear k r par_tmp CI_tmp
[par(1, :), CI(1, :, 1), CI(1, :, 2); par(2, :), CI(2, :, 1), CI(2, :, 2); par(3, :), CI(3, :, 1), CI(3, :, 2)]

% -340:
%    0.0039   12.1990  143.6397   -0.0123  -19.0271  136.5930    0.0201   43.4252  150.6865
%    0.0249   28.2646  146.2733   -0.0186   -4.0873  138.6268    0.0684   60.6165  153.9198
%    0.0109   15.5800  110.5960   -0.0182  -10.6801  104.9908    0.0401   41.8402  116.2012
% -370:
%    0.0055   16.5141  144.1475   -0.0189  -28.0735  136.2152    0.0299   61.1018  152.0798
%    0.0181   25.7481  145.3466   -0.0177   -7.2227  138.4723    0.0539   58.7189  152.2209
%    0.0130   17.0841  110.5151   -0.0196  -13.5493  105.1463    0.0456   47.7174  115.8839


% z_boundary = -340;
% % Pre_Process(z_boundary);
% load(['Presynaptic_coordinates_processed_', num2str(-z_boundary), '.mat'],...
    % 'Pre_L4E', 'Pre_L23E', 'Pre_L23I');

% DoubleGau_2 = @(r, kappa1e4, sigma_n10, sigma_b)...
    % ((kappa1e4 / 1e4) / (2 * pi * (sigma_n10 / 10)^2)) * exp(- r.^2 / (2 * (sigma_n10 / 10)^2)) +...
    % ((1 - (kappa1e4 / 1e4)) / (2 * pi * sigma_b^2)) * exp(- r.^2 / (2 * sigma_b^2));
% DoubleGau_1 = @(r, kappa1e4, sigma_n10, sigma_b)...
    % DoubleGau_2(r, kappa1e4, sigma_n10, sigma_b) .* (2 * pi * r);
% %
% opt = statset('Display', 'off', 'MaxIter', 1800, 'MaxFunEvals', 1800,...
    % 'TolX', 1e-8, 'TolFun', 1e-8, 'TolBnd', 1e-8);
% % par_ic = [0.005 * 1e4, 10 * 10, 150;...
% %     0.01 * 1e4, 30 * 10, 150;...
% %     0.01 * 1e4, 20 * 10, 100];
% par_ic = [0.05 * 1e4, 20 * 10, 150;...
    % 0.05 * 1e4, 20 * 10, 150;...
    % 0.05 * 1e4, 20 * 10, 100];
% %
% namelist = {'L4E', 'L23E', 'L23I'}; N_type = length(namelist);
% par = zeros(N_type, 3); CI = zeros(N_type, 3, 2); %single_sigma = zeros(1, N_type);
% for k = 1: N_type
    % eval(['r = sqrt(Pre_', namelist{k}, '(:, 1) .^ 2 + Pre_', namelist{k}, '(:, 2) .^ 2);']);
    % [par_tmp, CI_tmp] = mle(r, 'pdf',...
        % @(r, kappa1e4, sigma_n10, sigma_b)DoubleGau_1(r, kappa1e4, sigma_n10, sigma_b),...
        % 'start', par_ic(k, :), 'LowerBound', [0 0 0], 'UpperBound', [5000, 5000, 5000], 'options', opt);
    % par(k, :) = par_tmp; CI(k, :, 1) = CI_tmp(1, :); CI(k, :, 2) = CI_tmp(2, :);
    % %single_sigma(k) = sqrt(mean(r .^ 2) / 2);
% end
% clear k r par_tmp CI_tmp
% par(:, 1) = par(:, 1) / 1e4; par(:, 2) = par(:, 2) / 10;
% CI(:, 1, :) = CI(:, 1, :) / 1e4; CI(:, 2, :) = CI(:, 2, :) / 10; 

% [par(1, :), CI(1, :, 1), CI(1, :, 2);...
 % par(2, :), CI(2, :, 1), CI(2, :, 2);...
 % par(3, :), CI(3, :, 1), CI(3, :, 2)]




%% Figure
r1 = -500: 1: 500; r2 = 0: 1: 500;
%
Clb = [0.5 0 0.5; 1 0 0; 0 0 1];
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.8, 0.75]);
%
subplot(1, 2, 1); hold on;
l = zeros(1, 6);
for k = 1: 3
    l(2 * k - 1) = plot(r1, DoubleGau_2(r1, par(k, 1), par(k, 2), par(k, 3)), 'Color', Clb(k, :), 'LineWidth', 1.5);
end
for k = 1: 3
    l(2 * k) = plot(r1, -ones(size(r1)), 'Color', [1 1 1]);
end
legend(l, {'L2/3 Exc. <- L4 Exc.', '(feedforward exc.)',...
    'L2/3 Exc. <- L2/3 Exc.' ,'(recurrent exc.)',...
    'L2/3 Exc. <- L2/3 Inh.' ,'(recurrent inh.)'}, 'Fontsize', 12);
grid on;
axis([-500 500 0 2.1e-5]); set(gca, 'XTick', -500: 50: 500, 'YTick', 0: 0.25e-5: 2e-5);
set(gca, 'XTickLabel', {'-500', '', '-400', '', '-300', '', '-200', '', '-100', '', '0', '', '100', '', '200', '', '300', '', '400', '', '500'});
xlabel('Horizontal radial distance (\mum)'); ylabel('(p.d.f.)');
title(['(L2/3 and L4 boundary: ', num2str(z_boundary), ' \mum)'], 'FontWeight', 'normal');
%
subplot(1, 2, 2); hold on;
l = zeros(1, 6);
for k = 1: 3
    l(2 * k - 1) = plot(r2, DoubleGau_1(r2, par(k, 1), par(k, 2), par(k, 3)), 'Color', Clb(k, :), 'LineWidth', 1.5);
end
for k = 1: 3
    l(2 * k) = plot(r2, -ones(size(r2)), 'Color', [1 1 1]);
end
legend(l, {'L2/3 Exc. <- L4 Exc.', ['\kappa = ', num2str(par(1, 1), '%.3f'),...
    ', \sigmaN = ', num2str(par(1, 2), '%.3f'), ', \sigmaB = ', num2str(par(1, 3), '%.2f'), ' (\mum)'],...
    'L2/3 Exc. <- L2/3 Exc.', ['\kappa = ', num2str(par(2, 1), '%.3f'),...
    ', \sigmaN = ', num2str(par(2, 2), '%.3f'), ', \sigmaB = ', num2str(par(2, 3), '%.2f'), ' (\mum)'],...
    'L2/3 Exc. <- L2/3 Inh.', ['\kappa = ', num2str(par(3, 1), '%.3f'),...
    ', \sigmaN = ', num2str(par(3, 2), '%.3f'), ', \sigmaB = ', num2str(par(3, 3), '%.2f'), ' (\mum)']}, 'Fontsize', 12);
grid on;
axis([0 500 0 7e-3]); set(gca, 'XTick', 0: 50: 500, 'YTick', 0: 0.5e-3: 7e-3);
xlabel('Horizontal radial distance (\mum)'); ylabel('(p.d.f.)');
%
pause(1);
print(gcf, '-dpng', ['Conn_fit_', num2str(-z_boundary), '.png']);
close;


