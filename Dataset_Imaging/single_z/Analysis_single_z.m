%% This is for .mat file of individual z and individual choice of neuropil factor.
%% Attention: All pairs of OS ROIs will be preserved for calculation of correlation.
  % In some of non-hand-drawn-ROI datasets, for neuron pairs with
  % tuning correlation > 0.95 (there are cases of 0.99!), I manually removed one of them.
input_filename = 'result.mat';
output_filename = 'results_analyzed.mat';

%% Basic parameters for Y18
N_frame = [24, 6, 18];
pixel2um = 1;    % 1 um / 1 pixel
frame_rate = 2.18;    % Hz (frame / s)
N_trial = 10;
theta_stim = linspace(0, 330, 12); N_theta = length(theta_stim);
Thr_Pvalue_anova = 0.05;
Thr_Responsive = 0.1;
%
% Choice of distance coordinate for correlation ~ dist plot.
d_BinStep = 5;    % interval of distance bin, 5 um.
N_Bin = 201;
d_BinCenter = [([0: 1: (N_Bin - 2)]' + 0.5) * d_BinStep; 1250];
% For d_BinStep = 5 --   2.5, 7.5, 12.5, ...
d_BinEdge = [[0: 1: (N_Bin - 1)] * d_BinStep, 1500];
% For d_BinStep = 5 -- 0,   5,   10,   15, ...
% But we'll drop the 1st one (from 0 to 5, center at 2.5), since there will be few pairs.

%% Load
load(input_filename,...
    'cellNum', 'Intensity_raw', 'Intensity_neuropil', 'order', 'neuropilFactor', 'xy', 'stdImg');
N_neuron = cellNum;
ROI_size = size(stdImg) * pixel2um;    % in um.

%% Position of each ROI
ROI_ij = cell(1, N_neuron);    % (i, j) of all pixels for each ROI, in um.
for neuron_k = 1: N_neuron, ROI_ij{neuron_k} = xy{neuron_k}(:, [2, 1]) * pixel2um; end
% 1st column, x, actually j; 2nd column, y, actually i (not N - y).
ROI_centeroid = NaN(N_neuron, 2);    % (i, j) of center of mass for each ROI, in um.
for neuron_k = 1: N_neuron, ROI_centeroid(neuron_k, :) = mean(ROI_ij{neuron_k}, 1); end

%% Calculate dFF0.
% Intensity_raw and Intensity_neuropil are all what we need.
% baseline_neuropil: the average value of the first 150 min of Intensity_neuropil.
% For any alpha you want, the value "F" is
     % Intensity = Intensity_raw - alpha * bsxfun(@minus, Intensity_neuropil, baseline_neuropil);
% And "F0", baseline, make sure after got F:
    % histogram into 50 slices. Find the average value of the slice with most numbers.
% Then dFF0 = (F - F0) / F0.
baseline_neuropil = NaN(1, N_neuron);
for neuron_k = 1: N_neuron
    neuropil_tmp = sort(Intensity_neuropil(:, neuron_k));
    baseline_neuropil(neuron_k) = mean(neuropil_tmp(1: 150));
end
dF_neuropil = bsxfun(@minus, Intensity_neuropil, baseline_neuropil);
Intensity = Intensity_raw - bsxfun(@times, neuropilFactor, dF_neuropil);    % F
baseline = NaN(1, N_neuron);    % F0
for neuron_k = 1: N_neuron
    Intensity_k = Intensity(:, neuron_k);
    [Intensity_hist, Intensity_BinCtr] = hist(Intensity_k, 50);
    Intensity_dBin = Intensity_BinCtr(2) - Intensity_BinCtr(1);
    f0_MaxBin = find(Intensity_hist == max(Intensity_hist)); f0_MaxBin = f0_MaxBin(1);
    % The frame number of each neuron's most frequent intensity (idx).
    f0_group = find((Intensity_k > Intensity_BinCtr(f0_MaxBin) - Intensity_dBin / 2) &...
        (Intensity_k < Intensity_BinCtr(f0_MaxBin) + Intensity_dBin / 2));
    baseline(neuron_k) = mean(Intensity_k(f0_group), 'omitnan');
end
dF = bsxfun(@minus, Intensity, baseline);
dFF0 = bsxfun(@rdivide, dF, baseline);

%% Turn raw dFF0 into tuning for each trial.
% order(:, 1: 2): Each 120 (N_theta * N_trial) rows is one trial. 1st column is the direction used, 2nd is trial index.
% order(:, 3: 4): sortrows(order(:, 1: 2), 1). Use anyone of the two.
dFF0_Trials = NaN(N_neuron, N_theta, N_trial);
Trial_count = zeros(N_theta, 1);
for trial_k = 1: (N_theta * N_trial)
    idx1 = (trial_k - 1) * N_frame(1) + 1; theta_k = order(idx1, 1);
    theta_idx = (theta_k / (theta_stim(2) - theta_stim(1))) + 1;
    frame_idx = order((idx1 + N_frame(2)): (idx1 + N_frame(2) + N_frame(3) - 1), 2);
    Trial_count(theta_idx) = Trial_count(theta_idx) + 1;
    dFF0_Trials(:, theta_idx, Trial_count(theta_idx)) = mean(dFF0(frame_idx, :)', 2);
end
tmp = isnan(dFF0_Trials);
if (~isequal(Trial_count, N_trial * ones(N_theta, 1))) | (sum(tmp(:))), error('Trial count error.'); end
clear tmp
%
% Trial average and standard error
dFF0_mean = mean(dFF0_Trials, 3);
dFF0_se = std(dFF0_Trials, [], 3) / sqrt(N_trial); 

%% One-way ANOVA analysis: the p-value of flatten tuning curve.
Pvalue_anova = NaN(N_neuron, 1);
for neuron_k = 1: N_neuron
    Pvalue_anova(neuron_k) = anova1(squeeze(dFF0_Trials(neuron_k, :, :))',...
        {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l'}, 'off');
end

%% Fit preferred direction
Tuning_par_fit = NaN(N_neuron, 9);
% 9 elements are: theta_pref, sigma_pref, sigma_oppo, R_pref, R_oppo, R_offset, SSE, RSQ, exitflag
theta_stim_ext = [theta_stim, 360]; theta_stim_ext_interp = 0: 7.5: 360;
Tuning_curve_func = @(theta_stim, theta_pref, sigma_pref, sigma_oppo, R_pref, R_oppo, R_offset)...
    R_pref * exp(- dtheta_dir(theta_stim, theta_pref) .^ 2 / (2 * sigma_pref^2)) +...
    R_oppo * exp(- dtheta_dir(theta_stim, (theta_pref - pi)) .^ 2 / (2 * sigma_oppo^2)) + R_offset;
options = optimoptions('lsqnonlin', 'Display', 'none', 'MaxFunEvals', 1200, 'MaxIter', 1200);
for neuron_k = 1: N_neuron
    Tuning_curve = [dFF0_mean(neuron_k, :), dFF0_mean(neuron_k, 1)];
    [R_pref_0, idx] = max(Tuning_curve); theta_pref_0 = theta_stim_ext(idx);
    Tuning_curve = interp1(theta_stim_ext, Tuning_curve, theta_stim_ext_interp, 'spline');
    %
    Err = @(par) Tuning_curve_func(theta_stim_ext_interp * pi/180,...
        par(1), par(2), par(3), par(4), par(5), par(6)) - Tuning_curve;
    IniVal = [theta_pref_0, pi/4, pi/4, R_pref_0, 0, 0];
    par_lb = [- pi, 0, 0, 0, 0, -2 * R_pref_0];
    par_ub = [3 * pi, pi/2, pi/2, 2 * R_pref_0, 2 * R_pref_0, 2 * R_pref_0];
    [par, ~, ~, Tuning_par_fit(neuron_k, 9), ~, ~, ~] = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
    if par(4) < par(5)    % exchange pref and oppo
        tmp = par(4); par(4) = par(5); par(5) = tmp;
        tmp = par(2); par(2) = par(3); par(3) = tmp;
        if par(1) > pi, par(1) = par(1) - pi; else, par(1) = par(1) + pi; end; clear tmp
    end
    if par(1) < 0, par(1) = par(1) + 2 * pi; end
    if par(1) > 2 * pi, par(1) = par(1) - 2 * pi; end
    %
    Tuning_par_fit(neuron_k, 1: 6) = par;
    Tuning_par_fit(neuron_k, 7) = sum((Tuning_curve_func(theta_stim * pi/180,...
        par(1), par(2), par(3), par(4), par(5), par(6)) - dFF0_mean(neuron_k, :)) .^ 2);
        Tuning_par_fit(neuron_k, 8) = 1 - Tuning_par_fit(neuron_k, 7) /...
            (sum((dFF0_mean(neuron_k, :) - mean(dFF0_mean(neuron_k, :))) .^ 2) * ((N_theta - 1) / N_theta));
end
%
Pref_direction = Tuning_par_fit(:, 1) * 180/pi;    % in degree rather than rad
Pref_orientation = mod(Pref_direction, 180);
isGoodFitting = ((Tuning_par_fit(:, 7) < 0.4) & (Tuning_par_fit(:, 8) > 0.6));

%% If neurons are OS
isNonFlatten = (Pvalue_anova < Thr_Pvalue_anova);
isResponsive = (max(dFF0_mean, [], 2) > Thr_Responsive);
isOS = NaN(N_neuron, 2);
isOS(:, 1) = isNonFlatten & isResponsive;
isOS(:, 2) = isNonFlatten & isResponsive & isGoodFitting;
% When calculating tuning curve correlation ~ distance, I did not include isGoodFitting,
  % or there will be much fewer pairs for nearby distance.
% But when calculating \Delta\theta_pref ~ distance, we have to include isGoodFitting.

%% Tuning correlation / \Delta\theta_pref ~ distance
% For tuning correlation
ROI_centeroid_OS_1 = ROI_centeroid(isOS(:, 1) == 1, :);
dFF0_mean_OS = dFF0_mean(isOS(:, 1) == 1, :);
[SigCorr_1d, dist_1d_sigcorr] = SignalCorr_Dist(dFF0_mean_OS, ROI_centeroid_OS_1);
% Each element of represents one neuron pair.
  % For multiple z and multiple datasets, simply concatenate them.
  % But do not concatenate dFF0_mean_OS and ROI_centeroid_OS,
  % or you'll mix neuron pairs at different z / different datasets!
%
% For \Delta\theta_pref
ROI_centeroid_OS_2 = ROI_centeroid(isOS(:, 2) == 1, :);
Pref_orientation_OS = Pref_orientation(isOS(:, 2) == 1);
[dtheta_pref_1d, dist_1d_dthetapref] = dtheta_pref_Dist(Pref_orientation_OS, ROI_centeroid_OS_2);
%
[SigCorr_mean, SigCorr_se, Num_Bin_SigCorr] = histogram_mean_sem(SigCorr_1d, dist_1d_sigcorr, d_BinEdge);
[dtheta_pref_mean, dtheta_pref_se, Num_Bin_dtheta_pref] = histogram_mean_sem(dtheta_pref_1d, dist_1d_dthetapref, d_BinEdge);
%
idx_valid = find(Num_Bin_SigCorr >= 2);    % We'll drop intervals with 0 or only 1 single pair.
d_BinCenter_SigCorr = d_BinCenter(idx_valid);
SigCorr_mean = SigCorr_mean(idx_valid);
SigCorr_se = SigCorr_se(idx_valid);
Num_Bin_SigCorr = Num_Bin_SigCorr(idx_valid);
%
idx_valid = find(Num_Bin_dtheta_pref >= 2);    % We'll drop intervals with 0 or only 1 single pair.
d_BinCenter_dtheta_pref = d_BinCenter(idx_valid);
dtheta_pref_mean = dtheta_pref_mean(idx_valid);
dtheta_pref_se = dtheta_pref_se(idx_valid);


save(output_filename,...
    'N_neuron', 'ROI_size', 'ROI_ij', 'ROI_centeroid',...
    'dFF0_mean', 'dFF0_se', 'Tuning_curve_func', 'Tuning_par_fit',...
    'Pref_direction', 'Pref_orientation', 'isOS',...
    'd_BinCenter_SigCorr', 'SigCorr_mean', 'SigCorr_se', 'Num_Bin_SigCorr',...
    'd_BinCenter_dtheta_pref', 'dtheta_pref_mean', 'dtheta_pref_se', 'Num_Bin_dtheta_pref');


%% Plot
lwdth = 1.5; mksize = 12.5; cpsize = 7.5; txtsz = 15;
xylim = [0 100 -0.5 1.01]; xtickc = [0: 5: 100]; ytickc = -0.5: 0.1: 1;
%
figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
%
subplot(1, 2, 1); hold on;
errorbar(d_BinCenter_SigCorr, SigCorr_mean, SigCorr_se,...
    'Color', [1 0 0], 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
plot([0, 1000], [0 0], 'k--');
axis([0 100 -0.5 1]); set(gca, 'XTick', [0: 5: 100], 'YTick', [-0.5: 0.1: 1]);
axis square; grid on;
xlabel('Horizontal Cortical Distance (\mum)');
title('Tuning curve correlation', 'FontWeight', 'normal');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');
%
subplot(1, 2, 2); hold on;
errorbar(d_BinCenter_dtheta_pref, dtheta_pref_mean, dtheta_pref_se,...
    'Color', [1 0 0], 'LineWidth', lwdth, 'Marker', '.', 'MarkerSize', mksize, 'CapSize', cpsize);
plot([0, 1000], [45 45], 'k--');
axis([0 100 0 90]); set(gca, 'XTick', [0: 5: 100], 'YTick', [0: 5: 90]);
axis square; grid on; 
xlabel('Horizontal Cortical Distance (\mum)');
title('\Delta\theta', 'FontWeight', 'normal');
set(gca, 'FontSize', txtsz, 'box', 'off', 'TickDir', 'out');


