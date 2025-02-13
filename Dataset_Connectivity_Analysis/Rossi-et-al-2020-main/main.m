% The naming convention of the files matches the one for the files on Github (datasets 1-17). You have 2 types of .mat files:
% -VisualEX contains the data for the excitatory presynaptic neurons
% -VisualIN contains the data for the inhibitory neurons
% In each of the files, in addition to what's on github, you find the 2 variables you asked me for:
% - xyz: (dx, dy, z)  coordinate of each presynaptic neuron relative to the postsynaptic cell, in microns
% - deltaDir: difference in preferred direction of each presynaptic neuron with respect to the postsynaptic cell. The preferred direction is measured as the stimulus that drove the maximal response, hence the value is limited to the values deltaDir =  -180:30:150. 
% One important point: for now, there isn't an exact correspondence between the positions xyz in these files and the data you find on github in the spatialEX and spatialIN files on github. This is because of several reasons - these data come from the functional ROIs detected with Suite2P on longitudinal recordings, rather than from the structural analysis of 1 anatomical stack of images: not all neurons have visual responses, therefore some which are detected in the anatomical images, are not detected functionally; during longitudinal recordings, the brain shape or the recording pars may change, and some neurons may be duplicated or absent because they appear in different positions; some presynaptic neurons may appear and then disappear (toxicity of rabies); finally I am still working to try to align the two datasets as best as I can.
% With this in mind, if possible, I would use this data to test your assumption P(i, j) = P(r(i, j)) * P(Δθ_{pref.}(i, j)), but I would stick to the anatomical data for the previous analysis concerning the spatial probability of connections. 

% Random test 0.5 !
% Qs:
% More neurons in visual than spatial?
% What about IN?
%% Dependence on dr -- probably also dtheta !
%% > 0.5? > 0.05? ......
%% Other methods than p-values...


dir_0 = 'D:/Study/CompNeuro/Projects/1_Spatial_Clustering/Rossi-et-al-2020-main/Independence_data';
cd(dir_0);
%
Pre_TotE = [];
name_list = dir('visualEX*.mat'); name_list = {name_list.name};
for k = 1: length(name_list)
    load([dir_0, '/', name_list{k}], 'xyz', 'deltaDir');    % (Δx, Δy, -z)    % ranging from -180: 30: 150 deg
    r_k = sqrt(xyz(:, 1) .^ 2 + xyz(:, 2) .^ 2);
    Pre_TotE = [Pre_TotE; [r_k, deltaDir, xyz]];    % (r, dDir, Δx, Δy, -z)
end
clear xyz deltaDir k r_k name_list
cd ..
%
z_boundary = -370;
z_tmp = -Pre_TotE(:, 5);
idx_L4E = find((z_tmp < z_boundary) & (z_tmp >= -470));
idx_L23E = find((z_tmp <= -100) & (z_tmp >= z_boundary));
Pre_L23E = Pre_TotE(idx_L23E, :);
Pre_L4E = Pre_TotE(idx_L4E, :);


% Try again -- Figure of (Δx, Δy) scatters
dir1 = -180; dir2 = 150; ddir = 30;
dir_ctr = dir1: ddir: dir2; dir_edge = (dir1 - ddir / 2): ddir: (dir2 + ddir / 2); Ndir = length(dir_ctr);
%
%
%
%
clb_hsv = hsv(Ndir) * 0.875;
color_L23E = NaN(size(Pre_L23E, 1), 3); color_L4E = NaN(size(Pre_L4E, 1), 3);
for k = 1: size(Pre_L23E, 1);
    color_L23E(k, :) = clb_hsv(round((Pre_L23E(k, 2) - dir1) / ddir + 1), :);
end
for k = 1: size(Pre_L4E, 1);
    color_L4E(k, :) = clb_hsv(round((Pre_L4E(k, 2) - dir1) / ddir + 1), :);
end
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 0.8]);
subplot(1, 2, 1); scatter(Pre_L23E(:, 3), Pre_L23E(:, 4), 20, color_L23E, 'filled');
title(['L2/3 Exc. <- L2/3 Exc. (n = ', num2str(size(Pre_L23E, 1)), ')'], 'FontWeight', 'normal');
subplot(1, 2, 2); scatter(Pre_L4E(:, 3), Pre_L4E(:, 4), 20, color_L4E, 'filled');
title(['L2/3 Exc. <- L4 Exc. (n = ', num2str(size(Pre_L4E, 1)), ')'], 'FontWeight', 'normal');
colormap(clb_hsv); clb = colorbar;
set(clb, 'Position', [0.92 0.11 0.02 0.815], 'YTick', [0.5: Ndir - 0.5] / Ndir, 'YTickLabel', dir_ctr);
clbtitle = get(clb, 'Title'); set(clbtitle, 'String', '\DeltaDir (deg.)');
for k = 1: 2
    subplot(1, 2, k); axis square; grid on; axis([-400 400 -400 400]);
    set(gca, 'XTick', -400: 100: 400, 'YTick', -400: 100: 400);
    xlabel('\Deltax (\mum)'); ylabel('\Deltay (\mum)');
end
pause(2); print(gcf, '-dpng', 'Scatter_L23E_L4E.png'); close;


dr_list = [10, 15, 20, 24, 30, 40, 48];
for rk = 1: length(dr_list)
% Chi-square test of independence
% hist2
dr = dr_list(rk); r1 = dr; r2 = 360;    % ignore 1st interval 0 ~ dr
r_ctr = r1: dr: r2; r_edge = (r1 - dr / 2): dr: (r2 + dr / 2); Nr = length(r_ctr);
%dir1 = -180; dir2 = 150; ddir = 30;
%dir_ctr = dir1: ddir: dir2; dir_edge = (dir1 - ddir / 2): ddir: (dir2 + ddir / 2); Ndir = length(dir_ctr);
%
p2_L23E = histcounts2(Pre_L23E(:, 1), Pre_L23E(:, 2), r_edge, dir_edge);    % (r, dir)
p2_L4E = histcounts2(Pre_L4E(:, 1), Pre_L4E(:, 2), r_edge, dir_edge);    % (r, dir)
%
% dir -> theta
theta1 = -60; theta2 = 90; dtheta = 30;
theta_ctr = theta1: dtheta: theta2; Ntheta = length(theta_ctr);
p2_L23E = p2_L23E(:, 1: 6) + p2_L23E(:, 7: 12);    % 0: 30: 150
p2_L23E = p2_L23E(:, [5 6 1 2 3 4]);    % -60: 30: 90
p2_L4E = p2_L4E(:, 1: 6) + p2_L4E(:, 7: 12);
p2_L4E = p2_L4E(:, [5 6 1 2 3 4]);
p2_L23E_ext = [p2_L23E(:, end) p2_L23E];
p2_L4E_ext = [p2_L4E(:, end) p2_L4E];
p2_L23E_ext_rownorm = bsxfun(@rdivide, p2_L23E_ext, sum(p2_L23E_ext, 2));
p2_L4E_ext_rownorm = bsxfun(@rdivide, p2_L4E_ext, sum(p2_L4E_ext, 2));
%
% Get the p-value of H0: independent
p_value_L23E = Chi_square_test(p2_L23E);
p_value_L4E = Chi_square_test(p2_L4E);


figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);  % 0.6
% 2D histogram
subplot(2, 4, 1); imagesc(p2_L23E_ext_rownorm);
title({'2D histogram of L2/3 Exc. input number count (each row norm. by max.)',...
    ['p-value of Chi-square test of independence: ', num2str(p_value_L23E, '%.3f')]}, 'FontWeight', 'normal');
subplot(2, 4, 3); imagesc(p2_L4E_ext_rownorm);
title({'2D histogram of L4 Exc. input number count (each row norm. by max.)',...
    ['p-value of Chi-square test of independence: ', num2str(p_value_L4E, '%.3f')]}, 'FontWeight', 'normal');
clb1 = colorbar; set(clb1, 'Position', [0.7 0.582 0.01 0.34]);
if rk == 7, set(clb1, 'Position', [0.7 0.602 0.01 0.305]); end
for k = [1 3]
subplot(2, 4, k); axis equal; axis([0.5, Ntheta + 1.5, 0.5, Nr + 0.5]);
set(gca, 'XTick', 1: Ntheta + 1, 'YTick', 1: Nr, 'XTickLabel', [-90 theta_ctr], 'YTickLabel', r_ctr);
xlabel('\Delta Pref. Orientation (deg.)'); ylabel('Horizontal distance (\mum)');
end
% 2D histogram
subplot(2, 4, 5); imagesc(p2_L23E_ext);
title('(Raw number count)', 'FontWeight', 'normal');
subplot(2, 4, 7); imagesc(p2_L4E_ext);
title('(Raw number count)', 'FontWeight', 'normal');
ymax = max([max(p2_L23E(:)), max(p2_L4E(:))]);
clb2 = colorbar; set(clb2, 'Position', [0.7 0.11 0.01 0.34], 'YTick', 0: ymax);
for k = [5 7]
subplot(2, 4, k); axis equal; axis([0.5, Ntheta + 1.5, 0.5, Nr + 0.5]);
set(gca, 'XTick', 1: Ntheta + 1, 'YTick', 1: Nr, 'XTickLabel', [-90 theta_ctr], 'YTickLabel', r_ctr);
xlabel('\Delta Pref. Orientation (deg.)'); ylabel('Horizontal distance (\mum)');
end
% Marginal, L2/3
subplot(4, 4, 2);
ptheta_L23E = sum(p2_L23E, 1);
plot([-90 theta_ctr], [ptheta_L23E(end) ptheta_L23E], 'marker', '.', 'markersize', 15);
axis([-90 90 0 120]); set(gca, 'XTick', -90: 30: 90, 'YTick', 0: 30: 120);
xlabel('\Delta Pref. Orientation (deg.)'); ylabel('Number count');
title('Marginal distribution of L2/3E \Delta\theta', 'FontWeight', 'normal');
subplot(4, 4, 6);
pr_L23E = sum(p2_L23E, 2)'; plot(r_ctr, pr_L23E, 'marker', '.', 'markersize', 15);
axis([0 360 0 120]); set(gca, 'XTick', 0: 60: 360, 'YTick', 0: 30: 120);
xlabel('Horizontal distance (\mum)'); ylabel('Number count');
title('Marginal distribution of L2/3E r', 'FontWeight', 'normal');
% Marginal, L4
subplot(4, 4, 4);
ptheta_L4E = sum(p2_L4E, 1);
plot([-90 theta_ctr], [ptheta_L4E(end) ptheta_L4E], 'marker', '.', 'markersize', 15);
axis([-90 90 0 120]); set(gca, 'XTick', -90: 30: 90, 'YTick', 0: 30: 120);
xlabel('\Delta Pref. Orientation (deg.)'); ylabel('Number count');
title('Marginal distribution of L4E \Delta\theta', 'FontWeight', 'normal');
subplot(4, 4, 8);
pr_L4E = sum(p2_L4E, 2)'; plot(r_ctr, pr_L4E, 'marker', '.', 'markersize', 15);
axis([0 360 0 120]); set(gca, 'XTick', 0: 60: 360, 'YTick', 0: 30: 120);
xlabel('Horizontal distance (\mum)'); ylabel('Number count');
title('Marginal distribution of L4E r', 'FontWeight', 'normal');
%
pause(2); print(gcf, '-dpng', ['P_value_L23E_L4E_dr', num2str(dr), '.png']); close;
%
end




% Random test......
N_rand = 5000;
N_L23E = size(Pre_L23E, 1); N_L4E = size(Pre_L4E, 1);
p_value_L23E_rand = NaN(1, N_rand); p_value_L4E_rand = NaN(1, N_rand);
for k = 1: N_rand
    p2_L23E_randk = histcounts2(Pre_L23E(randperm(N_L23E), 1),...
        Pre_L23E(randperm(N_L23E), 2), r_edge, dir_edge);
    p_value_L23E_rand(k) = Chi_square_test(p2_L23E_randk);
    p2_L4E_randk = histcounts2(Pre_L4E(randperm(N_L4E), 1),...
        Pre_L4E(randperm(N_L4E), 2), r_edge, dir_edge);
    p_value_L4E_rand(k) = Chi_square_test(p2_L4E_randk);
    %fprintf([num2str(k), '/1000.\n']);
end

[mean(p_value_L23E_rand), std(p_value_L23E_rand); 0.5, 1/sqrt(12)]
[mean(p_value_L4E_rand), std(p_value_L4E_rand); 0.5, 1/sqrt(12)]

% UNIFORM DISTRIBUTION
