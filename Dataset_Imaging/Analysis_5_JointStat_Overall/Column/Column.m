if ispc, dirD = 'D:'; elseif isunix, dirD = '/media/DATA1'; end
addpath(genpath([dirD, '/Study/CompNeuro/Projects/Functions_simul/']));
dir0 = [dirD, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging'];
dir_local = [dir0, '/Analysis_5_JointStat_Overall/Column'];
load([dir_local, '/General_information_Z.mat']);


for dr_threshold = [5, 7.5, 10]
%dr_threshold = 7.5;

dz_center_1_individual = cell(1, Dataset_N_tot);
SigCorr_mean_individual = cell(1, Dataset_N_tot);
SigCorr_se_individual = cell(1, Dataset_N_tot);
dz_center_1_overall = cell(1, length(Dataset_type));
SigCorr_mean_overall= cell(1, length(Dataset_type));
SigCorr_se_overall = cell(1, length(Dataset_type));
%
dz_center_2_individual = cell(1, Dataset_N_tot);
dtheta_pref_mean_individual = cell(1, Dataset_N_tot);
dtheta_pref_se_individual = cell(1, Dataset_N_tot);
dz_center_2_overall = cell(1, length(Dataset_type));
dtheta_pref_mean_overall= cell(1, length(Dataset_type));
dtheta_pref_se_overall = cell(1, length(Dataset_type));


for i = 1: length(Dataset_type)
SigCorr_column_tot = []; dz_column_1_tot = [];
dtheta_pref_column_tot = []; dz_column_2_tot = [];
%
for dataset_k = (Dataset_idx_edge(i) + 1): Dataset_idx_edge(i + 1)
z_value_list = Dataset_z{dataset_k};
factor_k = neuropilFactor_best_idx(dataset_k);
Coor3_OS_1_tot = []; dFF0_mean_OS_tot = [];
Coor3_OS_2_tot = []; Pref_orientation_OS_tot = [];
%
for z_k = 1: length(z_value_list)
    z = z_value_list(z_k);
    savename = [Dataset_type{i}, '_', Dataset_name{dataset_k}, '_z', num2str(z)];
    load([dir0, '/Analysis_1_Individual_z/', savename, '_s2p_NpSize30_Ana1.mat'],...
        'N_neuron', 'ROI_centeroid', 'dFF0_mean', 'Pref_orientation', 'isOS');
    % remove "invalid" neurons
    idx_remove_k = idx_remove{dataset_k}{z_k};    % Comes from neuropil factor auto. No difference...
    idx_valid = 1: N_neuron; idx_valid(idx_remove_k) = NaN; idx_valid = idx_valid(~isnan(idx_valid));
    ROI_centeroid = ROI_centeroid(idx_valid, :);
    dFF0_mean = dFF0_mean(idx_valid, :, factor_k);
    Pref_orientation = Pref_orientation(idx_valid, factor_k) * (pi / 180);
    isOS = isOS(idx_valid, :, factor_k);
    % pick up OS
    idxOS_1 = find(isOS(:, 1) == 1); idxOS_2 = find(isOS(:, 2) == 1);
    ROI_centeroid_OS_1 = ROI_centeroid(idxOS_1, :); dFF0_mean_OS = dFF0_mean(idxOS_1, :);
    ROI_centeroid_OS_2 = ROI_centeroid(idxOS_2, :); Pref_orientation_OS = Pref_orientation(idxOS_2);
    Coor3_OS_1 = [ROI_centeroid_OS_1, z * ones(size(ROI_centeroid_OS_1, 1), 1)];    % (i, j, z) in micron
    Coor3_OS_2 = [ROI_centeroid_OS_2, z * ones(size(ROI_centeroid_OS_2, 1), 1)];
    % record
    Coor3_OS_1_tot = [Coor3_OS_1_tot; Coor3_OS_1];
    dFF0_mean_OS_tot = [dFF0_mean_OS_tot; dFF0_mean_OS];
    Coor3_OS_2_tot = [Coor3_OS_2_tot; Coor3_OS_2];
    Pref_orientation_OS_tot = [Pref_orientation_OS_tot; Pref_orientation_OS];
end
clear idx_remove_k idx_valid ROI_centeroid dFF0_mean Pref_orientation isOS idxOS_1 idxOS_2
clear ROI_centeroid_OS_1 ROI_centeroid_OS_2 Coor3_OS_1 dFF0_mean_OS Coor3_OS_2 Pref_orientation_OS
%
%% SigCorr, each dataset
di = bsxfun(@minus, Coor3_OS_1_tot(:, 1), Coor3_OS_1_tot(:, 1)');
dj = bsxfun(@minus, Coor3_OS_1_tot(:, 2), Coor3_OS_1_tot(:, 2)');
dr = sqrt(di .^ 2 + dj .^ 2); dr = triu_new(dr, 0, 1); clear di dj
dz = abs(bsxfun(@minus, Coor3_OS_1_tot(:, 3), Coor3_OS_1_tot(:, 3)')); dz = triu_new(dz, 0, 1);
SigCorr = corrcoef(dFF0_mean_OS_tot'); SigCorr = triu_new(SigCorr, 0, 1);
idx_column = find(dr <= dr_threshold);
SigCorr_column = SigCorr(idx_column); dz_column = dz(idx_column);
clear idx_column SigCorr dr dz Coor3_OS_1_tot dFF0_mean_OS_tot
%
dz_center_k = sort(unique(dz_column)); Nz = length(dz_center_k);
SigCorr_mean = NaN(1, Nz); SigCorr_se = NaN(1, Nz);
for m = 1: Nz
    SigCorr_m = SigCorr_column(abs(dz_column - dz_center_k(m)) < 0.1);
    SigCorr_mean(m) = mean(SigCorr_m);
    err = SigCorr_m - SigCorr_mean(m); SigCorr_se(m) = std(err) / sqrt(length(err));
end
clear m SigCorr_m err Nz
%
SigCorr_column_tot = [SigCorr_column_tot; SigCorr_column];
dz_column_1_tot = [dz_column_1_tot; dz_column];
dz_center_1_individual{dataset_k} = dz_center_k;
SigCorr_mean_individual{dataset_k} = SigCorr_mean;
SigCorr_se_individual{dataset_k} = SigCorr_se;
clear SigCorr_column dz_column dz_center_k SigCorr_mean SigCorr_se
%
%% dtheta_pref, each dataset
di = bsxfun(@minus, Coor3_OS_2_tot(:, 1), Coor3_OS_2_tot(:, 1)');
dj = bsxfun(@minus, Coor3_OS_2_tot(:, 2), Coor3_OS_2_tot(:, 2)');
dr = sqrt(di .^ 2 + dj .^ 2); dr = triu_new(dr, 0, 1); clear di dj
dz = abs(bsxfun(@minus, Coor3_OS_2_tot(:, 3), Coor3_OS_2_tot(:, 3)')); dz = triu_new(dz, 0, 1);
dtheta_pref = abs(bsxfun(@dtheta, Pref_orientation_OS_tot, Pref_orientation_OS_tot')) * (180 / pi);
dtheta_pref = triu_new(dtheta_pref, 0, 1);
idx_column = find(dr <= dr_threshold);
dtheta_pref_column = dtheta_pref(idx_column); dz_column = dz(idx_column);
clear idx_column dtheta_pref dr dz Coor3_OS_2_tot Pref_orientation_OS_tot
%
dz_center_k = sort(unique(dz_column)); Nz = length(dz_center_k);
dtheta_pref_mean = NaN(1, Nz); dtheta_pref_se = NaN(1, Nz);
for m = 1: Nz
    dtheta_pref_m = dtheta_pref_column(abs(dz_column - dz_center_k(m)) < 0.1);
    dtheta_pref_mean(m) = mean(dtheta_pref_m);
    err = dtheta_pref_m - dtheta_pref_mean(m); dtheta_pref_se(m) = std(err) / sqrt(length(err));
end
clear m dtheta_pref_m err Nz
%
dtheta_pref_column_tot = [dtheta_pref_column_tot; dtheta_pref_column];
dz_column_2_tot = [dz_column_2_tot; dz_column];
dz_center_2_individual{dataset_k} = dz_center_k;
dtheta_pref_mean_individual{dataset_k} = dtheta_pref_mean;
dtheta_pref_se_individual{dataset_k} = dtheta_pref_se;
clear dtheta_pref_column dz_column dz_center_k dtheta_pref_mean dtheta_pref_se
%
end
%
%% SigCorr, overall per dataset
dz_center_k = sort(unique(dz_column_1_tot)); Nz = length(dz_center_k);
SigCorr_mean = NaN(1, Nz); SigCorr_se = NaN(1, Nz);
for m = 1: Nz
    SigCorr_m = SigCorr_column_tot(abs(dz_column_1_tot - dz_center_k(m)) < 0.1);
    SigCorr_mean(m) = mean(SigCorr_m);
    err = SigCorr_m - SigCorr_mean(m); SigCorr_se(m) = std(err) / sqrt(length(err));
end
dz_center_1_overall{i} = dz_center_k';
SigCorr_mean_overall{i} = SigCorr_mean;
SigCorr_se_overall{i} = SigCorr_se;
clear m SigCorr_m err Nz SigCorr_column_tot dz_column_1_tot dz_center_k SigCorr_mean SigCorr_se
%
%% dtheta pref, overall per dataset
dz_center_k = sort(unique(dz_column_2_tot)); Nz = length(dz_center_k);
dtheta_pref_mean = NaN(1, Nz); dtheta_pref_se = NaN(1, Nz);
for m = 1: Nz
    dtheta_pref_m = dtheta_pref_column_tot(abs(dz_column_2_tot - dz_center_k(m)) < 0.1);
    dtheta_pref_mean(m) = mean(dtheta_pref_m);
    err = dtheta_pref_m - dtheta_pref_mean(m); dtheta_pref_se(m) = std(err) / sqrt(length(err));
end
dz_center_2_overall{i} = dz_center_k;
dtheta_pref_mean_overall{i} = dtheta_pref_mean;
dtheta_pref_se_overall{i} = dtheta_pref_se;
clear m dtheta_pref_m err Nz dtheta_pref_column_tot dz_column_1_tot
clear dz_center_k dtheta_pref_mean dtheta_pref_se
%
end





Exp2 = @(x, A, lambda, b) A * exp(- (x / (2 * lambda)).^2) + b;
x2 = linspace(0, 100, 1001);
IniVal = [1, 20, 0]; par_lb = [0, 0, 0]; par_ub = [1000, 1000, 1000];
options = optimoptions('lsqnonlin', 'Display', 'none', 'MaxFunEvals', 1200, 'MaxIter', 1200);
%
lwdth = 1.5; mksize = 12.5; cpsize = 7.5; txtsz = 15;
clb = [1 0.5 0; 1 0 0; 0 0 1; 0.5 0 1];
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.55, 1]); hold on;
suptitle({['Signal Correlation ~ Vertical Cortical Distance (\Deltaz)'],...
        ['(Horizontal radius = ', num2str(dr_threshold), ' \mum)']}, 8, 0.95);
l = zeros(1, length(Dataset_type));
lgdtxt_0 = {'L2/3 cytosolic GCaMP', 'L2/3 nuclear GCaMP', 'L4 cytosolic GCaMP', 'L4 soma GCaMP'};
lgdtxt = cell(1, length(Dataset_type));
for i = 1: length(Dataset_type)
    x_data = dz_center_1_overall{i};
    y_data = SigCorr_mean_overall{i};
    y_data_se = SigCorr_se_overall{i};
    l(i) = errorbar(x_data, y_data, y_data_se,...
        'Color', clb(i, :), 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
    %
    if dr_threshold == 7.5
        Err = @(par) (Exp2(x_data, par(1), par(2), par(3)) - y_data) ./ y_data_se;
        [fitpar, ~, residual, exitflag, ~, ~, Jacobian] = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
        CI = nlparci(fitpar, residual, 'jacobian', Jacobian); fitpar_err = ((CI(:, 2) - CI(:, 1)) / 2)';
        plot(x2, Exp2(x2, fitpar(1), fitpar(2), fitpar(3)), 'LineStyle', '-.', 'Color', clb(i, :));
        lgdtxt{i} = [lgdtxt_0{i}, ' (\lambda = ', num2str(fitpar(2), '%.2f'),...
            ' \pm ', num2str(fitpar_err(2), '%.2f'), ' \mum)'];
    end
end
plot([0, 101], [0 0], 'k--');
axis([0, 101, -0.1, 1.01]); set(gca, 'XTick', 0: 10: 100, 'YTick', [-0.1, 0: 0.1: 1]);
% set(gca, 'XTick', 0: 20: 100, 'YTick', [-0.1, 0: 0.25: 1]);
axis square; grid on; 
legend(l, lgdtxt, 'FontSize', txtsz);
xlabel('Vertical Cortical Distance, \Deltaz (\mum)');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
pause(2); print(gcf, '-dpng', [dir_local, '/Figures/SigCorr_overall_', num2str(dr_threshold, '%.1f'), '.png']);
if dr_threshold == 7.5
    print(gcf, '-depsc2', '-r1500', [dir_local, '/Figures/SigCorr_overall_', num2str(dr_threshold, '%.1f'), '.eps']);
end
close;

figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.55, 1]); hold on;
suptitle({['\Delta\theta\_pref. (deg.) ~ Vertical Cortical Distance (\Deltaz)'],...
        ['(Horizontal radius = ', num2str(dr_threshold), ' \mum)']}, 8, 0.95);
l = zeros(1, length(Dataset_type));
for i = 1: length(Dataset_type)
    l(i) = errorbar(dz_center_2_overall{i}, dtheta_pref_mean_overall{i}, dtheta_pref_se_overall{i},...
        'Color', clb(i, :), 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
end
plot([0, 101], [45 45], 'k--');
axis([0, 101, 0, 61]); set(gca, 'XTick', 0: 10: 100, 'YTick', 0: 5: 60);
% set(gca, 'XTick', 0: 20: 100, 'YTick', [0 20 40 45]);
axis square; grid on; 
lgdtxt = {'L2/3 cytosolic GCaMP', 'L2/3 nuclear GCaMP', 'L4 cytosolic GCaMP', 'L4 soma GCaMP'};
legend(l, lgdtxt, 'FontSize', txtsz, 'Location', 'northwest');
xlabel('Vertical Cortical Distance, \Deltaz (\mum)');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
pause(2); print(gcf, '-dpng', [dir_local, '/Figures/dtheta_pref_overall_', num2str(dr_threshold, '%.1f'), '.png']);
if dr_threshold == 7.5
    print(gcf, '-depsc2', '-r1500', [dir_local, '/Figures/dtheta_pref_overall_', num2str(dr_threshold, '%.1f'), '.eps']);
end
close;


end


