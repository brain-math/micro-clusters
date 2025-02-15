Parset_i = 2;


load([pwd, '/Parameters/JN_', num2str(Parset_i), '.mat'], 'JN_ef', 'JN_ee', 'JN_ei');
load([pwd, '/Results/Results_Total_', num2str(Parset_i), '.mat'], 'Pref_theta_F', 'd_tot',...
    'Nf_sqrt', 'Ne_sqrt', 'Ni_sqrt', 'Pref_theta_E', 'Pref_theta_I', 'sigma_micron', 'kappa');
d_f = d_tot / Nf_sqrt; d_e = d_tot / Ne_sqrt; d_i = d_tot / Ni_sqrt;
N_post_sample = [4000, 4000, 2500];


dist_ctr = cell(1, 3); JN_d_1d_pdf = cell(1, 3); JN_d_par = cell(1, 3);
dtheta_ctr = cell(1, 3); JN_dtheta_pdf = cell(1, 3); JN_dtheta_par = cell(1, 3);
%
[dist_ctr{1}, JN_d_1d_pdf{1}, JN_d_par{1}, dtheta_ctr{1}, JN_dtheta_pdf{1}, JN_dtheta_par{1}] =...
    Wiring_analysis(JN_ef, N_post_sample(1), d_f, Pref_theta_E, Pref_theta_F);
clear JN_ef
[dist_ctr{2}, JN_d_1d_pdf{2}, JN_d_par{2}, dtheta_ctr{2}, JN_dtheta_pdf{2}, JN_dtheta_par{2}] =...
    Wiring_analysis(JN_ee, N_post_sample(2), d_e, Pref_theta_E, Pref_theta_E);
clear JN_ee
[dist_ctr{3}, JN_d_1d_pdf{3}, JN_d_par{3}, dtheta_ctr{3}, JN_dtheta_pdf{3}, JN_dtheta_par{3}] =...
    Wiring_analysis(JN_ei, N_post_sample(3), d_i, Pref_theta_E, Pref_theta_I);
clear JN_ei

Gau2 = @(r, sigma) exp(- r .^ 2 / (2 * sigma ^ 2)) / (2 * pi * sigma ^ 2);
DoubleGau_2 = @(r, A, kappa, sigma_n, sigma_b)...
    A * ( kappa * Gau2(r, sigma_n) + (1 - kappa) * Gau2(r, sigma_b) );
DoubleGau_1 = @(r, A, kappa, sigma_n, sigma_b)...
    DoubleGau_2(r, A, kappa, sigma_n, sigma_b) .* (2 * pi * r);
Gau0 = @(dtheta, A_dtheta, sigma_dtheta, a_dtheta)...
    A_dtheta * exp(- dtheta.^2 / (2 * sigma_dtheta ^ 2)) + a_dtheta;
dist_0 = 0.1: 0.1: 350; theta_0 = 0: 0.1: 90;
%
JN_d_1d_pdf_fit = cell(1, 3); JN_d_1d_pdf_before = cell(1, 3);
JN_dtheta_1 = cell(1, 3); JN_dtheta_1_par = cell(1, 3);
JN_dtheta_pdf_fit = cell(1, 3); JN_dtheta_1_fit = cell(1, 3);
for k = 1: 3
    JN_d_1d_pdf_fit{k} = DoubleGau_1(dist_0, JN_d_par{k}(1),...
        JN_d_par{k}(2), JN_d_par{k}(3), JN_d_par{k}(4));
    JN_d_1d_pdf_before{k} = DoubleGau_1(dist_0, 1,...
        kappa(k), sigma_micron(2, k), sigma_micron(1, k));
    %
    A0 = (JN_dtheta_par{k}(1) + JN_dtheta_par{k}(3));
    JN_dtheta_1{k} = JN_dtheta_pdf{k} / A0;
    JN_dtheta_1_par{k} = [JN_dtheta_par{k}(2), JN_dtheta_par{k}(3) / A0];
    JN_dtheta_pdf_fit{k} = Gau0(theta_0, JN_dtheta_par{k}(1),...
        JN_dtheta_par{k}(2), JN_dtheta_par{k}(3));
    JN_dtheta_1_fit{k} = JN_dtheta_pdf_fit{k} / A0;
end
%
save([pwd, '/Figures/Wiring_analysis_', num2str(Parset_i), '.mat']);






figure; set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.82, 1]);
namelist_title = {'EF', 'EE', 'EI'};
%
for k = 1: 3
%
l = zeros(1, 2);
subplot(2, 3, k); hold on;
plot(dist_ctr{k}, JN_d_1d_pdf{k}, 'Marker', '.', 'Color', 'r', 'MarkerSize', 10, 'LineStyle', 'none');
l(1) = plot(dist_0, JN_d_1d_pdf_before{k}, 'b--');
l(2) = plot(dist_0, JN_d_1d_pdf_fit{k}, 'r');
legend(l, {['Before: (', num2str(kappa(k), '%.3f'), ', ', num2str(sigma_micron(2, k), '%.3f'), ', ',...
    num2str(sigma_micron(1, k), '%.1f'), ')'], ['After: (', num2str(JN_d_par{k}(2), '%.3f'), ', ',...
    num2str(JN_d_par{k}(3), '%.3f'), ', ', num2str(JN_d_par{k}(4), '%.1f'), ')']});
axis([0 300 0 0.0125]); set(gca, 'XTick', 0: 50: 300, 'YTick', 0: 0.0025: 0.0125);
axis square; grid on; xlabel('Distance (\mum)');
title(['J_{', namelist_title{k}, '}(d) (normalized p.d.f.)'], 'fontweight', 'normal', 'fontsize', 12);
%
subplot(2, 3, k + 3); hold on;
plot(dtheta_ctr{k}, JN_dtheta_1{k}, 'b');
plot(theta_0, JN_dtheta_1_fit{k}, 'r--');
legend(['\sigma(\Delta\theta) = ', num2str(JN_dtheta_1_par{k}(1), '%.2f'),...
    '^0, a = ', num2str(JN_dtheta_1_par{k}(2), '%.2f')]);
axis([0 90 0 1.2]); set(gca, 'XTick', 0: 15: 90, 'YTick', 0: 0.2: 1.2);
axis square; grid on; xlabel('|\Delta\theta| (deg.)');
title(['J_{', namelist_title{k}, '}(\Delta\theta) (normalized p.d.f.)'], 'fontweight', 'normal', 'fontsize', 12);
%
end
%
pause(2); print(gcf, '-dpng', [pwd, '/Figures/Wiring_analysis_', num2str(Parset_i), '.png']); close;

