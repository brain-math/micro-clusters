Parset_i = 2;
rng('shuffle');

load([pwd, '/Results/Results_Total_', num2str(Parset_i), '.mat'],...
    'L23_FR', 'L23_Isyn', 'Ne');
FR_E = L23_FR(1: Ne, :); FR_I = L23_FR(Ne + 1: end, :);
IsynTot_E = L23_Isyn(1: Ne, :, end); IsynTot_I = L23_Isyn(Ne + 1: end, :, end);
%
N_sample = 750; N_total = [length(FR_E(:)), length(FR_I(:))];
idx_sample = NaN(2, N_sample);
for k = 1: 2, tmp = randperm(N_total(k)); idx_sample(k, :) = sort(tmp(1: N_sample)); end
%
IsynTot_E_scatter = IsynTot_E(idx_sample(1, :)); FR_E_scatter = FR_E(idx_sample(1, :));
IsynTot_I_scatter = IsynTot_I(idx_sample(2, :)); FR_I_scatter = FR_I(idx_sample(2, :));
clear L23_FR L23_Isyn N_total k tmp IsynTot_E IsynTot_I
%
SSN = @(a, n, thr, I) (a * (I - thr) .^ n) .* ((sign(I - thr) + 1) / 2);
E_err = @(par) SSN(par(1), par(2), par(3), IsynTot_E_scatter) - FR_E_scatter;
I_err = @(par) SSN(par(1), par(2), par(3), IsynTot_I_scatter) - FR_I_scatter;
IniVal = [35, 1, -0.3]; par_lb = [0, 0, -1000]; par_ub = [1000, 1000, 1000];
options = optimoptions('lsqnonlin', 'Display', 'none', 'MaxFunEvals', 1200, 'MaxIter', 1200);
[ParE, ~, ~, ~] = lsqnonlin(E_err, IniVal, par_lb, par_ub, options); aE = ParE(1); nE = ParE(2); thrE = ParE(3);
[ParI, ~, ~, ~] = lsqnonlin(I_err, IniVal, par_lb, par_ub, options); aI = ParI(1); nI = ParI(2); thrI = ParI(3);
clear E_err I_err IniVal par_lb par_ub options ParE ParI
%
OPyE = mean(FR_E(:)); OPxE = (OPyE / aE) ^ (1 / nE) + thrE;
Ge = aE * nE * ((OPyE / aE) ^ ((nE - 1) / nE));
OPyI = mean(FR_I(:)); OPxI = (OPyI / aI) ^ (1 / nI) + thrI;
Gi = aI * nI * ((OPyI / aI) ^ ((nI - 1) / nI));
clear FR_E FR_I


lwdth = 2; ftsize = 12; titleftsize = 14;
%
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
Isyn = linspace(-0.25, 0.75, 1001); ymax = 25; dy = 2.5;
%
subplot(1, 2, 1); hold on; grid on;
scatter(IsynTot_E_scatter, FR_E_scatter, 10, 'fill',...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'c'); 
plot(Isyn, SSN(aE, nE, thrE, Isyn), 'Color', 'r', 'LineWidth', lwdth);
scatter(OPxE, OPyE, 50, 'fill', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
set(gca, 'YTick', [0: dy: ymax]); axis([min(Isyn) max(Isyn) 0 ymax]);
xlabel('Input (V/s)', 'FontSize', titleftsize); ylabel('Firing Rate (Hz)', 'FontSize', titleftsize);
title({['f (Hz) = ', num2str(aE, '%.2f'), '\cdot[I (mV/ms) + ', num2str(- thrE, '%.2f'),...
    ']_+^{', num2str(nE, '%.2f'), '}'], ['ge = ', num2str(Ge, '%.3f'), ' (V^{-1})']},...
    'FontSize', titleftsize, 'FontWeight', 'normal');
%
subplot(1, 2, 2); hold on; grid on;
scatter(IsynTot_I_scatter, FR_I_scatter, 10, 'fill',...
    'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'c');
plot(Isyn, SSN(aI, nI, thrI, Isyn), 'Color', 'b', 'LineWidth', lwdth);
scatter(OPxI, OPyI, 50, 'fill', 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
set(gca, 'YTick', [0: dy: ymax]); axis([min(Isyn) max(Isyn) 0 ymax]);
xlabel('Input (V/s)', 'FontSize', titleftsize); ylabel('Firing Rate (Hz)', 'FontSize', titleftsize);
title({['f (Hz) = ', num2str(aI, '%.2f'), '\cdot[I (mV/ms) + ', num2str(- thrI, '%.2f'),...
    ']_+^{', num2str(nI, '%.2f'), '}'], ['gi = ', num2str(Gi, '%.3f'), ' (V^{-1})']},...
    'FontSize', titleftsize, 'FontWeight', 'normal');
%
pause(2); print(gcf, '-dpng', [pwd, '/Figures/Gain_', num2str(Parset_i), '.png']); close;

