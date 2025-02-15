if ispc, dirD = 'D:'; elseif isunix, dirD = '/media/DATA1'; end
addpath(genpath([dirD, '/Study/CompNeuro/Projects/Functions_simul/']));
dir0 = [dirD, '/Study/CompNeuro/Projects/Micro-clustering/Dataset_Imaging'];
load([dir0, '/General_information.mat']);
%
theta_stim_ext = [theta_stim, 360]; theta_stim_ext_interp = 0: 7.5: 360;
Tuning_curve_func = @(theta_stim, theta_pref, sigma_pref, sigma_oppo, R_pref, R_oppo, R_offset)...
     R_pref * exp(- dtheta_dir(theta_stim, theta_pref) .^ 2 / (2 * sigma_pref^2)) +...
     R_oppo * exp(- dtheta_dir(theta_stim, (theta_pref - pi)) .^ 2 / (2 * sigma_oppo^2)) + R_offset;
options = optimoptions('lsqnonlin', 'Display', 'none', 'MaxFunEvals', 1200, 'MaxIter', 1200);


for i = 1: length(Dataset_type)
for dataset_k = (Dataset_idx_edge(i) + 1): Dataset_idx_edge(i + 1)
z_value_list = Dataset_z{dataset_k};
for z_k = 1: length(z_value_list)
    raw_data_name = [Dataset_type{i}, '_', Dataset_name{dataset_k}, '_z',...
        num2str(z_value_list(z_k)), '_s2p_NpSize30.mat'];
    load([dir0, '/RawData/', raw_data_name], 'cellNum', 'Intensity_raw', 'Intensity_neuropil',...
        'order', 'neuropilFactor_auto', 'ROI_ij', 'stdImg', 'gOSI', 'xls_info');
    N_frame_k = N_frame{dataset_k}(z_k, :);
    %
    % Each z (i.e. each file in RawData) has its individual result .mat file. Following terms will be saved:
        % N_neuron: scalar
        % ROI_size: (row, column) in micron.
        % ROI_centeroid: (N_neuron, 2). Each row: (i, j) in micron. Conversion to (x, y): x = j; y = ROI_size(1) - y;
        % dFF0_mean: (N_neuron, N_theta, N_neuropilFactor)
        % dFF0_se: (N_neuron, N_theta, N_neuropilFactor)
        % Tuning_curve_func: function
        % Tuning_par_fit: (N_neuron, 9, N_neuropilFactor):
            % theta_pref, sigma_pref, sigma_oppo, R_pref, R_oppo, R_offset, SSE, RSQ, exitflag
        % Pref_direction, Pref_orientation: (N_neuron, N_neuropilFactor)
        % isOS: (N_neuron, 2, N_neuropilFactor). Standard 1 or 2.
            % 1 for correlation analysis, 2 for pref. orientation analysis.
        % neuropilFactor_auto_avg
    %
    N_neuron = cellNum; clear cellNum
    %
    pixel2um_k = pixel2um(dataset_k); ROI_size = size(stdImg) * pixel2um_k;
    %
    ROI_centeroid = NaN(N_neuron, 2);    % From pixel to um.
    for neuron_k = 1: N_neuron
        ROI_centeroid(neuron_k, :) = mean(ROI_ij{neuron_k}, 1) * pixel2um_k;
    end
    %
    dFF0_mean = NaN(N_neuron, N_theta, N_neuropilFactor);
    dFF0_se = NaN(N_neuron, N_theta, N_neuropilFactor);
    Tuning_par_fit = NaN(N_neuron, 9, N_neuropilFactor);
    Pref_direction = NaN(N_neuron, N_neuropilFactor);
    Pref_orientation = NaN(N_neuron, N_neuropilFactor);
    isOS = NaN(N_neuron, 2, N_neuropilFactor);
    %
    % Now we'll calculate dFF0.
    % Intensity_raw and Intensity_neuropil are the same. That's all what we need.
    % baseline_neuropil: the average value of the first 150 min of Intensity_neuropil.
    % For any alpha you want, the value "F" is
        % Intensity = Intensity_raw - alpha * bsxfun(@minus, Intensity_neuropil, baseline_neuropil);
    % And "F0", baseline, make sure after got F:
        % histogram into 50 slices. Find the average value of the slice with most numbers.
    % So dFF0 = (F - F0) / F0.
    baseline_neuropil = NaN(1, N_neuron);
    for neuron_k = 1: N_neuron
        neuropil_tmp = sort(Intensity_neuropil(:, neuron_k));
        baseline_neuropil(neuron_k) = mean(neuropil_tmp(1: 150));
    end
    Intensity_neuropil_2 = bsxfun(@minus, Intensity_neuropil, baseline_neuropil);
    %
    for factor_k = 1: N_neuropilFactor
        if factor_k <= length(neuropilFactor_manual)
            Intensity = Intensity_raw - neuropilFactor_manual(factor_k) * Intensity_neuropil_2;
        else
            Intensity = Intensity_raw - bsxfun(@times, neuropilFactor_auto, Intensity_neuropil_2);
        end
        %
        baseline = NaN(1, N_neuron);
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
        %
        dF = bsxfun(@minus, Intensity, baseline); dFF0 = bsxfun(@rdivide, dF, baseline);
        % Raw dFF0 -> tuning et each trial.
        % order(:, 1: 2): Each 120 (N_theta * N_trial) rows is one trial. 1st column is the direction used, 2nd is trial index.
        % order(:, 3: 4): sortrows(order(:, 1: 2), 1). Use anyone of the two.
        dFF0_Trials = NaN(N_neuron, N_theta, N_trial);
        Trial_count = zeros(N_theta, 1);
        for trial_k = 1: (N_theta * N_trial)
            idx1 = (trial_k - 1) * N_frame_k(1) + 1; theta_k = order(idx1, 1);
            theta_idx = (theta_k / (theta_stim(2) - theta_stim(1))) + 1;
            frame_idx = order((idx1 + N_frame_k(2)): (idx1 + N_frame_k(2) + N_frame_k(3) - 1), 2);
            Trial_count(theta_idx) = Trial_count(theta_idx) + 1;
            dFF0_Trials(:, theta_idx, Trial_count(theta_idx)) = mean(dFF0(frame_idx, :)', 2);
        end
        tmp = isnan(dFF0_Trials);
        if (~isequal(Trial_count, N_trial * ones(N_theta, 1))) | (sum(tmp(:)))
            error('Trial count error.');
        end
        clear tmp
        dFF0_mean_k = mean(dFF0_Trials, 3); dFF0_mean(:, :, factor_k) = dFF0_mean_k;
        dFF0_se(:, :, factor_k) = std(dFF0_Trials, [], 3) / sqrt(N_trial); 
        %
        % Fit preferred orientation
            % theta_pref, sigma_pref, sigma_oppo, R_pref, R_oppo, R_offset, SSE, RSQ, exitflag
        for neuron_k = 1: N_neuron
            Tuning_curve = [dFF0_mean_k(neuron_k, :), dFF0_mean_k(neuron_k, 1)];
            [R_pref_0, idx] = max(Tuning_curve); theta_pref_0 = theta_stim_ext(idx);
            Tuning_curve = interp1(theta_stim_ext, Tuning_curve, theta_stim_ext_interp, 'spline');
            %
            Err = @(par) Tuning_curve_func(theta_stim_ext_interp * pi/180,...
                par(1), par(2), par(3), par(4), par(5), par(6)) - Tuning_curve;
            IniVal = [theta_pref_0, pi/4, pi/4, R_pref_0, 0, 0];
            par_lb = [- pi, 0, 0, 0, 0, -2 * R_pref_0];
            par_ub = [3 * pi, pi/2, pi/2, 2 * R_pref_0, 2 * R_pref_0, 2 * R_pref_0];
            [par, ~, ~, Tuning_par_fit(neuron_k, 9, factor_k), ~, ~, ~] = lsqnonlin(Err, IniVal, par_lb, par_ub, options);
            par = real(par);    % JUST FOR L4_cytosolic_GCaMP_Y24_x0y0_z280, neuron_k = 130; Weird, but imag in 1e-6.
            if par(4) < par(5)
                tmp = par(4); par(4) = par(5); par(5) = tmp;
                tmp = par(2); par(2) = par(3); par(3) = tmp;
                if par(1) > pi, par(1) = par(1) - pi; else, par(1) = par(1) + pi; end
                clear tmp
            end
            if par(1) < 0, par(1) = par(1) + 2 * pi; end
            if par(1) > 2 * pi, par(1) = par(1) - 2 * pi; end
            %
            Tuning_par_fit(neuron_k, 1: 6, factor_k) = par;
            Tuning_par_fit(neuron_k, 7, factor_k) = sum((Tuning_curve_func(theta_stim * pi/180,...
                par(1), par(2), par(3), par(4), par(5), par(6)) - dFF0_mean_k(neuron_k, :)) .^ 2);
            Tuning_par_fit(neuron_k, 8, factor_k) = 1 - Tuning_par_fit(neuron_k, 7, factor_k) /...
                (sum((dFF0_mean_k(neuron_k, :) - mean(dFF0_mean_k(neuron_k, :))) .^ 2) * ((N_theta - 1) / N_theta));
        end
        Pref_direction(:, factor_k) = Tuning_par_fit(:, 1, factor_k) * 180/pi;
        Pref_orientation(:, factor_k) = mod(Pref_direction(:, factor_k), 180);
        isGoodFitting = ((Tuning_par_fit(:, 7, factor_k) < 0.4) & (Tuning_par_fit(:, 8, factor_k) > 0.6));
        %
        % One-way ANOVA analysis: the p-value of flatten tuning curve.
        Pvalue_anova = NaN(N_neuron, 1);
        for neuron_k = 1: N_neuron
            Pvalue_anova(neuron_k) = anova1(squeeze(dFF0_Trials(neuron_k, :, :))',...
                {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l'}, 'off');
        end
        %
        % If neurons are OS
        isNonFlatten = (Pvalue_anova < Thr_Pvalue_anova);
        isResponsive = (max(dFF0_mean_k, [], 2) > Thr_Responsive);
        isOS(:, 1, factor_k) = isNonFlatten & isResponsive;
        isOS(:, 2, factor_k) = isNonFlatten & isResponsive & isGoodFitting;
    end
    %
    savename = [Dataset_type{i}, '_', Dataset_name{dataset_k}, '_z',...
        num2str(z_value_list(z_k)), '_s2p_NpSize30_Ana1.mat'];
    save([dir0, '/Analysis_1_Individual_z/', savename], 'N_neuron', 'ROI_size', 'ROI_centeroid',...
        'dFF0_mean', 'dFF0_se', 'Tuning_curve_func', 'Tuning_par_fit', 'Pref_direction', 'Pref_orientation', 'isOS');
    fprintf([savename, ' saved.\n']);
end
end
end




