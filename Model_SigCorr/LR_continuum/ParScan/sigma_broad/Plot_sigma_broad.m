if ispc, dir_0 = 'D:'; elseif isunix, dir_0 = '/media/DATA1'; end
addpath(genpath([dir_0, '/Study/CompNeuro/Projects/Functions_simul/']));
dir_local = [dir_0, '/Study/CompNeuro/Projects/Micro-clustering/Model_SigCorr/LR_continuum'];
cd(dir_local); addpath([dir_local, '/ParScan']);
rng('shuffle');


scan_mode = 2;

load([dir_local, '/ParScan/sigma_broad/sigma_broad_mode', num2str(scan_mode), '.mat']);
%
switch scan_mode
case 1
y_min = -0.005; y_max = 0.025; x_max = 300;
y_txt = [0.14, 0.115, 0.09] * (y_max / 0.16);
dxy_full = 1; lambda_lb = 400; lambda_ub = 700; A_max = 12.5;
case 2
y_min = -0.02; y_max = 0.1; x_max = 300;
y_txt = [0.14, 0.115, 0.09] * (y_max / 0.16);
dxy_full = 1; lambda_lb = 400; lambda_ub = 700; A_max = 12.5;
end

figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 1]); hold on;
y_plot = 1: 2: Ny; x_plot = 1: 2: Nx;
Ny_plot = length(y_plot); Nx_plot = length(x_plot);
for i_plot = 1: Ny_plot
for j_plot = 1: Nx_plot
    i = y_plot(i_plot); j = x_plot(j_plot);
    k = (i_plot - 1) * Nx_plot + j_plot;
    subaxis(Ny_plot, Nx_plot, k, 'Spacing', 0.005, 'Padding', 0, 'Margin', 0.04);
    %
    plot(d_micron, squeeze(Cov_L23E(Ny + 1 - i, j, :))); hold on;
    plot([0 x_max], [0 0], 'k--');
    %
    axis square; axis([0 x_max y_min y_max]);
    set(gca, 'XTick', 0: (x_max / 6): x_max, 'YTick', [y_min: -y_min: y_max]);
    grid on; axftsize = 5; ax = gca; ax.FontSize = axftsize;
    if ~((i_plot == Ny_plot) & (j_plot == 1)), set(gca, 'XTickLabel', [], 'YTickLabel', []); end
    if j_plot == 1, if i_plot == round((Ny_plot + 1) / 2)
        ylabel({y_name, num2str(y_list(Ny + 1 - i))}, 'FontSize', 10);
    else
        ylabel(num2str(y_list(Ny + 1 - i)), 'FontSize', 10);
    end; end
    if i_plot == Ny_plot, if j_plot == round((Nx_plot + 1) / 2)
        xlabel({num2str(x_list(j)), x_name}, 'FontSize', 10);
    else
        xlabel(num2str(x_list(j)), 'FontSize', 10);
    end; end
end
end
%
pause(2); print(gcf, '-dpng', [dir_local,...
    '/ParScan/sigma_broad/Plot_curves_sigma_broad_mode', num2str(scan_mode), '.png']);
close;




A_raw = Cov_L23E_par(:, :, 2); lambda_raw = Cov_L23E_par(:, :, 1);
% [min(A_raw(:)), max(A_raw(:))]
% [min(lambda_raw(:)), max(lambda_raw(:))]
[X, Y] = meshgrid(x_list, y_list);
x_list_full = 0: dxy_full: max(x_list); y_list_full = 0: dxy_full: max(y_list);
[Xq, Yq] = meshgrid(x_list_full, y_list_full);
lambda = interp2(X, Y, lambda_raw, Xq, Yq);
A = interp2(X, Y, A_raw, Xq, Yq);
%
clb = jet(101) * 0.875;    %clb = clb(end: -1: 1, :);
lambda_value = (lambda - lambda_lb) / (lambda_ub - lambda_lb);
lambda_value(lambda_value < 0) = 0; lambda_value(lambda_value > 1) = 1;
color_idx = round(lambda_value * 100) + 1;
%
A_lb = 0; A_ub = A_max;
A_value = (A - A_lb) / (A_ub - A_lb);
A_value(A_value < 0) = 0; A_value(A_value > 1) = 1;
%
lambda_clr = ones(size(lambda, 1), size(lambda, 2), 3);    % white for default, NaN.
for i = 1: size(lambda, 1)
for j = 1: size(lambda, 2)
    try    % avoid NaN
    color_tmp = clb(color_idx(i, j), :);
    color_real = ones(1, 3) - A_value(i, j) * (ones(1, 3) - color_tmp);
    lambda_clr(i, j, :) = reshape(color_real, [1 1 3]);
    end
end
end
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 1]); hold on;
imagesc(x_list_full, y_list_full, lambda_clr); axis xy;
colormap(clb); clb = colorbar; axis square; caxis([lambda_lb lambda_ub]);
clb.Ticks = linspace(lambda_lb, lambda_ub, 7);    % [lambda_lb: 2.5: lambda_ub];
axis([x_list(1), x_list(end), y_list(1), y_list(end)]);
set(gca, 'XTick', x_list(1): 25: x_list(end), 'YTick', y_list(1): 25: y_list(end));
xlabel(x_name); ylabel(y_name);
axftsize = 15; ax = gca; ax.FontSize = axftsize;
set(gca, 'Layer', 'top');
%
pause(2); print(gcf, '-dpng', [dir_local,...
    '/ParScan/sigma_broad/Plot_colormap_sigma_broad_mode', num2str(scan_mode), '.png']);
close;


