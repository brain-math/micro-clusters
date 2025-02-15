Parset_i = 2;
load([pwd, '/Results/Results_Total_', num2str(Parset_i), '.mat'], 's_sample', 'Ne', 'Nrec');
load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i), '.mat'], 'param');
T_off = param.T_off; T_on = param.T_on;
T_label = [0, T_off, T_off + T_on];%, 2 * T_off + T_on, 2 * (T_off + T_on)];
clear param

figure;
mksize = 0.5;
scatter(s_sample(:, 1), s_sample(:, 2), mksize, 'k', 'filled'); axis ij;
axis([0 T_label(end) 1 Nrec]);
set(gca, 'XTick', T_label, 'YTick', [1 Ne Nrec]);
xlabel('Time (ms)'); ylabel('Neuron #');
%
pause(2); print(gcf, '-dpng', [pwd, '/Figures/Raster_eg_', num2str(Parset_i), '.png']); close;

