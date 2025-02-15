psp = @(t, tau_r, tau_d) (exp(- t / tau_r) - exp(- t / tau_d)) / (tau_r - tau_d);
t = linspace(0, 150, 1001);

epsp_fast = psp(t, 1, 5);
ipsp_fast = psp(t, 1, 8);
epsp_mixture = 0.7 * psp(t, 1, 5) + 0.3 * psp(t, 10, 100);
ipsp_mixture = 0.7 * psp(t, 1, 8) + 0.3 * psp(t, 10, 100);

figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 0.5]);
subplot(1, 2, 1); hold on;
plot(t, epsp_fast, 'r');
plot(t, epsp_mixture, 'm');
legend('Fast', '0.7 Fast + 0.3 Slow');
grid on;
title({'EPSP, \tau_r^{fast} = 1, \tau_d^{fast} = 10 (ms)',...
    '\tau_r^{slow} = 5, \tau_d^{slow} = 100 (ms)'}, 'fontweight', 'normal');
subplot(1, 2, 2); hold on;
plot(t, ipsp_fast, 'b');
plot(t, ipsp_mixture, 'color', [0 0.5 1]);
legend('Fast', '0.7 Fast + 0.3 Slow');
grid on;
title({'IPSP, \tau_r^{fast} = 1, \tau_d^{fast} = 10 (ms)',...
    '\tau_r^{slow} = 8, \tau_d^{slow} = 100 (ms)'}, 'fontweight', 'normal');
%
pause(2); print(gcf, '-dpng', 'psp_shape_0.3.png'); close;

