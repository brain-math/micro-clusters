d_BinCenter = [7.5: 5: 502.5]';
N_distBin = length(d_BinCenter);
%
NF_alpha = 0.75;
%
group_name_list = {'L23_evoked_NF075_dth0', 'L23_evoked_NF075_dth45',...
    'L23_evoked_NF075_dth90', 'L23_evoked_NF075_dth135'};%, 'L23_evoked_NF075_all'};
filename_list = {'L23cytosolic_pairCorr_fDist_dth0', 'L23cytosolic_pairCorr_fDist_dth45',...
    'L23cytosolic_pairCorr_fDist_dth90', 'L23cytosolic_pairCorr_fDist_dth135'};
N_group = length(group_name_list);
%
DualExp2 = @(x, A1, A2, lambda1, lambda2, b)...
    A1 * exp(- (x / lambda1).^2) + A2 * exp(- (x / lambda2).^2) + b;
IniVal = [0.5, 0.05, 10, 200, 0];
par_lb = [0, 0, 0, 0, -10]; par_ub = [100, 100, 10000, 10000, 10];
options = optimoptions('lsqnonlin', 'Display', 'none',...
    'MaxFunEvals', 1200, 'MaxIter', 1200, 'StepTolerance', 1e-8);
%
PairCorr_scatter = cell(1, N_group);
Corr_avg = NaN(N_distBin, N_group + 1);
Corr_se = NaN(N_distBin, N_group + 1);
Corr_std = NaN(N_distBin, N_group);
Num_pair = NaN(N_distBin, N_group);
FitPar = NaN(N_group + 1, 5);
FitPar_err = NaN(N_group + 1, 5);


tic;
for group_k = 1: N_group 

filename = filename_list{group_k};
[~, dataset_name] = xlsfinfo([filename, '.xlsx']);
idx_start = find(strcmp(dataset_name, 'Dataset0'));
dataset_name = dataset_name(idx_start: end); clear idx_start
%
% Raw scatters for each group.
N_dataset = length(dataset_name);
PairCorr = cell(1, N_dataset);
for dataset_k = 1: N_dataset
    tmp = xlsread([filename, '.xlsx'], dataset_name{dataset_k});
    tmp = tmp(2: end, 2: end); tmp = tmp(1: N_distBin, :);
    PairCorr{dataset_k} = tmp;
end
PairCorr_scatter{group_k} = PairCorr;
fprintf(['File ', num2str(group_k), ' / ', num2str(N_group),...
    ' has been read in ', num2str(toc / 60, '%.2f'), ' min.\n']);
clear filename dataset_name dataset_k tmp
%
% Average of scatters within each dist bin.
Corr_avg_datasets = NaN(N_distBin, N_dataset);
Num_pair_datasets = NaN(N_distBin, N_dataset);
for dataset_k = 1: N_dataset
    Corr_avg_datasets(:, dataset_k) = mean(PairCorr{dataset_k}, 2, 'omitnan');
    Num_pair_datasets(:, dataset_k) = sum(~isnan(PairCorr{dataset_k}), 2);
end
%
Num_pair_tot = sum(Num_pair_datasets, 2);
Corr_avg_tot = sum(...    % weighted average
    Corr_avg_datasets .*...
    bsxfun(@rdivide, Num_pair_datasets, Num_pair_tot),...    % weight
    2, 'omitnan');
clear dataset_k Corr_avg_datasets Num_pair_datasets
%
% Standard error.
Corr_MSE = zeros(N_distBin, 1);
for dataset_k = 1: N_dataset
    Corr_MSE = Corr_MSE + sum(bsxfun(@minus,...
        PairCorr{dataset_k}, Corr_avg_tot) .^ 2, 2, 'omitnan');
end
Corr_std_tot = sqrt(Corr_MSE) ./ sqrt(Num_pair_tot - 1);
Corr_se_tot = Corr_std_tot ./ sqrt(Num_pair_tot);
clear dataset_k Corr_MSE
clear PairCorr N_dataset
%
% Save
Corr_avg(:, group_k) = Corr_avg_tot;
Corr_se(:, group_k) = Corr_se_tot;
Corr_std(:, group_k) = Corr_std_tot; clear Corr_std_tot
Num_pair(:, group_k) = Num_pair_tot; clear Num_pair_tot
%
fprintf(['File ', num2str(group_k), ' / ', num2str(N_group),...
    ' done in ', num2str(toc / 60, '%.2f'), ' min.\n']);
%
end


% Integrate all
Corr_avg(:, N_group + 1) = sum(Corr_avg(:, 1: N_group) .* bsxfun(@rdivide, Num_pair, sum(Num_pair, 2)), 2);
%
Corr_MSE_tot = zeros(N_distBin, 1);
for group_k = 1: N_group 
for dataset_k = 1: length(PairCorr_scatter{group_k})
    Corr_MSE_tot = Corr_MSE_tot + sum(bsxfun(@minus,...
        PairCorr_scatter{group_k}{dataset_k}, Corr_avg(:, end)) .^ 2, 2, 'omitnan');
end
end
clear group_k dataset_k
Corr_std_tot = sqrt(Corr_MSE_tot) ./ sqrt(sum(Num_pair, 2) - 1);
Corr_se(:, N_group + 1) = Corr_std_tot ./ sqrt(sum(Num_pair, 2));
clear Corr_MSE_tot Corr_std_tot


% Fitting
for group_k = 1: N_group + 1
idx_valid = [1: 3, 6: length(d_BinCenter)];    % outlier 4 and 5?
Err = @(par) (DualExp2(d_BinCenter(idx_valid), par(1), par(2),...
    par(3), par(4), par(5)) - Corr_avg(idx_valid, group_k));% ./ Corr_se(idx_valid, group_k);
[par, ~, residual, ~, ~, ~, Jacobian] = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
CI = nlparci(par, residual, 'jacobian', Jacobian);
FitPar(group_k, :) = par;
FitPar_err(group_k, :) = ((CI(:, 2) - CI(:, 1)) / 2)';
end
clear group_k idx_valid Err par residual Jacobian CI IniVal options par_ub par_lb Corr_avg_tot Corr_se_tot


save([pwd, '/Results_Noise_correlation_new.mat'],...
    'group_name_list', 'N_group', 'NF_alpha', 'd_BinCenter',  'N_distBin',...
    'PairCorr_scatter', 'Corr_avg', 'Corr_se', 'Corr_std',...
    'Num_pair', 'DualExp2', 'FitPar', 'FitPar_err');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([pwd, '/Results_Noise_correlation_new.mat']);

group_idx = [1, 4, 5];
%lgdname = {'L23 (evoked rsp., \Delta\theta ~ 0)', 'L23 (evoked rsp., \Delta\theta ~ 45)',...
%    'L23 (evoked rsp., \Delta\theta ~ 90)', 'L23 (evoked rsp., \Delta\theta ~ 135)', 'L23 (evoked rsp, all)'};
lgdname = {'L23 (evoked rsp., similar)', 'L23 (evoked rsp., opposite)', 'L23 (evoked rsp, all)'};
clr = [0 0.5 1; 1 0.5 0; 0 0 0];%[1 0 0; 2/3 0 1/3; 1/3 0 2/3; 0 0 1; 0 0 0];
lwdth = 1.5; mksize = 5; cpsize = 7.5; txtsz = 12; lgdtxtsz = 11;
x2 = linspace(0, 50, 501);
ymax = [0.5, 0.5, 0.5]; dy = [0.05, 0.05, 0.05];
%
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
%
subplot(1, 2, 1); hold on;
l = zeros(1, length(group_idx));
for group_i = 1: length(group_idx)
    group_k = group_idx(group_i);
    l(group_i) = plot(d_BinCenter, Corr_avg(:, group_k),...
        'Color', clr(group_i, :), 'LineWidth', lwdth);
end
axis([0 500 0 ymax(1)]); set(gca, 'XTick', 0: 50: 500, 'YTick', 0: dy(1): ymax(1));
axis square; grid on;
legend(l, lgdname, 'FontSize', lgdtxtsz);
% legend(l, lgdname{group_idx}, 'FontSize', lgdtxtsz);
set(gca, 'FontSize', txtsz);
xlabel('Horizontal distance (\mum)');
title('Noise Correlation (NF = 0.75)', 'FontWeight', 'normal');
%
subplot(1, 2, 2); hold on;
l = zeros(1, length(group_idx));
for group_i = 1: length(group_idx)
    group_k = group_idx(group_i);
    l(group_i) = errorbar(d_BinCenter, Corr_avg(:, group_k), Corr_se(:, group_k),...
        'Color', clr(group_i, :), 'LineWidth', lwdth, 'Marker', '.',...
        'MarkerSize', mksize, 'CapSize', cpsize);
    par = FitPar(group_k, :);
    plot(x2, DualExp2(x2, par(1), par(2), par(3), par(4), par(5)),...
        'Color', clr(group_i, :), 'LineStyle', '--');
end
axis([0 50 0 ymax(2)]); set(gca, 'XTick', 0: 5: 50, 'YTick', 0: dy(2): ymax(2));
axis square; grid on;
legend(l, lgdname, 'FontSize', lgdtxtsz);
set(gca, 'FontSize', txtsz);
xlabel('Horizontal distance (\mum)');
title({'(Zoom in)', 'C = A_b exp(- (d / \lambda_b)^2) + A_n exp(- (d / \lambda_n)^2) + c_0'},...
    'FontWeight', 'normal');

%
pause(1); print(gcf, '-dpng', 'Results_Noise_correlation_evoked_new.png'); close;

