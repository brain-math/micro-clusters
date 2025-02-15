%par_idx = [3, 4, 5; 6, 7, 8; 9, 10, 11];
par_idx = [12, 13, 14; 15, 16, 17; 18, 19, 20];
% Broad: (100, 150, 150), (150, 150, 150), (150, 100, 100);
% kappa: (0.025, 0.025, 0.025), (0.125, 0.025, 0.025), (0.125, 0.125, 0.125);    % [144, 146, 110; 4, 10, 10]
% narrow: (10, 4, 4), (4, 4, 4), (4, 10, 10);    % [0.075, 0.125, 0.125]

clr = [1 0 0; 0.85 0.85 0; 0 0.75 0];
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
for i = 1: 3
for j = 1: 3
    load([pwd, '/Results/Results_Total_', num2str(par_idx(i, j)), '.mat'],...
        'd_BinCenter_E', 'SigCorr_mean_E', 'SigCorr_se_E', 'Pref_theta_E');
    subplot(3, 4, (i - 1) * 4 + 1); hold on;
    if i == 1, plot(d_BinCenter_E, SigCorr_mean_E, 'color', clr(j, :));
    % smooth(SigCorr_mean_E, 15)
    else, errorbar(d_BinCenter_E, SigCorr_mean_E, SigCorr_se_E, 'color', clr(j, :)); end
    subplot(3, 4, (i - 1) * 4 + (j + 1));
    imagesc(Pref_theta_E); colormap(hsv); colorbar;
    title('---------', 'color', clr(j, :));
end
end
subplot(3, 4, 1); axis([0 600 -0.05 0.15]); axis square; grid on;
set(gca, 'XTick', 0: 150: 600, 'YTick', -0.05: 0.05: 0.15);
ylabel('sigma broad'); legend('(100, 150, 150)', '(150, 150, 150)', '(150, 100, 100)');
%
subplot(3, 4, 5); axis([0 50 -0.1 0.6]); axis square; grid on;
set(gca, 'XTick', 0: 10: 50, 'YTick', -0.1: 0.1: 0.6);
ylabel('kappa'); legend('(.025, .025, .025)', '(.125, .025, .025)', '(.125, .125, .125)', 'fontsize', 7);
title('fixed sigma: [144, 146, 110; 4, 10, 10]', 'fontweight', 'normal');
%
subplot(3, 4, 9); axis([0 50 -0.1 0.6]); axis square; grid on;
set(gca, 'XTick', 0: 10: 50, 'YTick', -0.1: 0.1: 0.6);
ylabel('sigma narrow'); legend('(10, 4, 4)', '(4, 4, 4)', '(4, 10, 10)');
title('fixed kappa: [0.075, 0.125, 0.125]', 'fontweight', 'normal');

%pause(2); print(gcf, '-dpng', [pwd, '/Figures/SigCorr_ParScan_spatial.png']); close;
pause(2); print(gcf, '-dpng', [pwd, '/Figures/SigCorr_ParScan_dthetadep.png']); close;

