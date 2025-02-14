% Use Ian's fig.
% Fig105_all.fig is all 3 groups, using edge 7.5: 5: ..., center 10: 5: ...
% Fig121_v3.fig for all cells, using edge 15: 10: ..., center 20: 10: ...

dir_0 = 'D:\Study\CompNeuro\Projects\Micro-clustering\Dataset_EPerturb_OneCellATime\IanUpdate_December2024';

lwdth = 1; mksize = 30; cpsize = 10; txtsz = 15;


close all;
openfig([dir_0, '\Fig105_all.fig']);
h_all = findobj(gca, 'type', 'errorbar');
%
n = length(get(h_all(1), 'XData'));
x = NaN(3, n); y = NaN(3, n); y_se = NaN(3, n); clr = NaN(3, 3);
for k = 1: 3    % g, r, k
    x(k, :) = get(h_all(k), 'XData');
    y(k, :) = get(h_all(k), 'YData');
    y_se(k, :) = get(h_all(k), 'YNegativeDelta');
    clr(k, :) = get(h_all(k), 'Color');
end
%
close all;
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 0.67]);
for k = 1: 3
    if k <= 2, idx1 = 2; elseif k == 3, idx1 = 1; end
    subplot(1, 3, k); hold on;
    errorbar(x(k, idx1: end), y(k, idx1: end), y_se(k, idx1: end), 'color', clr(k, :),...
        'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
    plot([0 50], [0 0], 'k--');
    plot([15 15], [-0.1 0.5], 'k--');
    xlim([0 50]); ylim([-0.1 0.5]);
    axis square; set(gca, 'box', 'on');
end
%
pause(2); print(gcf, '-dpng', 'Fig_7D.png');
print(gcf, '-depsc2', '-r1500', 'Fig_7D.eps');
close all;
clear h_all n x y y_se clr k idx1



close all;
openfig([dir_0, '\Fig121_v3.fig']);
h = findobj(gca, 'type', 'errorbar');
%
x = get(h, 'XData');
y = get(h, 'YData');
y_se = get(h, 'YNegativeDelta');
% clr = get(h, 'Color');
%
close all;
figure; hold on;
errorbar(x, y, y_se, 'color', 'k',...
    'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
plot([0 75], [0 0], 'k--');
xlim([0 75]); ylim([-0.05 0.1]);
axis square; set(gca, 'box', 'on');
%
pause(2); print(gcf, '-dpng', 'Fig_7E.png');
print(gcf, '-depsc2', '-r1500', 'Fig_7E.eps');
close all;
clear h x y y_se k

