function Pre_Process(z_boundary)
% z_boundary = -370;

%The file Presynaptic_coordinates.mat contains 2 cell arrays, inCoords and pcCoords.
%Each of the 17 entries in these arrays correspond to a different tracing experiment,
%  and contains a nPresynapticNeurons x 3 matrix of spatial coordinates XYZ in microns:
%  coordinates X and Y represent the distance of the presynaptic neurons relative to
%  their postsynaptic neuron, while coordinate Z represents the depth of neurons
%  relative to the cortical surface.
%  inCoords contains the coordinates of inhibitory presynaptic neurons, while
%  pcCoords contains the coordinates of excitatory presynaptic neurons.
%  Presynaptic_coordinates.mat also includes a vector postsynZ, which contains
%  the Z coordinate of the postsynaptic neurons in each experiment (X and Y coordinates
%  of the postsynaptic neurons are 0 in this framework).
%As you are interested in the layer distribution of neurons, I should mention that
%  on average, in my preparation, I estimated the boundary between L2/3 and L4
%  to be at 340um depth, and the boundary between L4 and L5 to be at 470 um depth. 
%BTW, All postsynZ are within L2/3 depth.

load('Presynaptic_coordinates.mat', 'postsynZ', 'rotInCoord', 'rotPcCoord');
postsynZ = - postsynZ;
pcCoord = rotPcCoord; clear rotPcCoord    % Δx, Δy, z of presynaptic neurons
inCoord = rotInCoord; clear rotInCoord
%
N_post = length(postsynZ);    % 17 post neurons
N_pre_each = zeros(2, N_post);    % Number of E and I pre neurons (2 rows) for each post neuron (17 columns)
for k = 1: N_post
    N_pre_each(1, k) = size(pcCoord{k}, 1);
    N_pre_each(2, k) = size(inCoord{k}, 1);
end
N_pre_tot = sum(N_pre_each, 2)';    % Totally 957 E and 543 I.
%
% Put the (Δx, Δy, Δz, z) information of all presynaptic neurons (957 E and 543 I) together.
Pre_TotE = []; Pre_TotI = [];    % 957 or 543 rows; 4 columns: Δx, Δy, Δz, z
for k = 1: N_post
    Pre_E_tmp = []; Pre_I_tmp = [];
    Pre_E_tmp(:, [1 2 4]) = pcCoord{k}; Pre_E_tmp(:, 4) = - Pre_E_tmp(:, 4);
    Pre_E_tmp(:, 3) = Pre_E_tmp(:, 4) - postsynZ(k);
    Pre_TotE = cat(1, Pre_TotE, Pre_E_tmp);
    Pre_I_tmp(:, [1 2 4]) = inCoord{k}; Pre_I_tmp(:, 4) = - Pre_I_tmp(:, 4);
    Pre_I_tmp(:, 3) = Pre_I_tmp(:, 4) - postsynZ(k);
    Pre_TotI = cat(1, Pre_TotI, Pre_I_tmp);
end
clear k Pre_E_tmp Pre_I_tmp pcCoord inCoord
%
% Pick up L2/3 or L4 presynaptic neurons
idx_L4E = find((Pre_TotE(:, 4) < z_boundary) & (Pre_TotE(:, 4) >= -470));
idx_L23E = find((Pre_TotE(:, 4) <= -100) & (Pre_TotE(:, 4) >= z_boundary));
idx_L4I = find((Pre_TotI(:, 4) < z_boundary) & (Pre_TotI(:, 4) >= -470));
idx_L23I = find((Pre_TotI(:, 4) <= -100) & (Pre_TotI(:, 4) >= z_boundary));
%
% The (Δx, Δy, Δz) information of: L4E, L4I, L2/3E, L2/3I presynaptic neurons.
% And sqrt(Δx^2 + Δy^2) is the r of post-pre horizontal cortical distance.
Pre_L23E = Pre_TotE(idx_L23E, 1: 3); Pre_L23I = Pre_TotI(idx_L23I, 1: 3);    % Δx, Δy, Δz
Pre_L4E = Pre_TotE(idx_L4E, 1: 3); Pre_L4I = Pre_TotI(idx_L4I, 1: 3);
clear idx_L23E idx_L23I idx_L4E idx_L4I

save(['Presynaptic_coordinates_processed_', num2str(-z_boundary), '.mat']);



%%%%%%%%%%%% 3D plot of (Δx, Δy, Δz) of L4E, L4I, L2/3E, L2/3I presynaptic neurons.
Clb = [1 0 0; 0 0 1; 0.5 0 0.5; 1 0 1];
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.6, 1]);
subplot(2, 2, 1); scatter3(Pre_L4E(:, 1), Pre_L4E(:, 2), Pre_L4E(:, 3),...
    10, 'fill', 'MarkerEdgeColor', Clb(3, :), 'MarkerFaceColor', Clb(3, :));
title(['L2/3 Exc. <- L4 Exc. (n = ', num2str(length(Pre_L4E)), ')'], 'FontWeight', 'normal');
subplot(2, 2, 2); scatter3(Pre_L4I(:, 1), Pre_L4I(:, 2), Pre_L4I(:, 3),...
    10, 'fill', 'MarkerEdgeColor', Clb(4, :), 'MarkerFaceColor', Clb(4, :));
title(['L2/3 Exc. <- L4 Inh. (n = ', num2str(length(Pre_L4I)), ')'], 'FontWeight', 'normal');
subplot(2, 2, 3); scatter3(Pre_L23E(:, 1), Pre_L23E(:, 2), Pre_L23E(:, 3),...
    10, 'fill', 'MarkerEdgeColor', Clb(1, :), 'MarkerFaceColor', Clb(1, :));
title(['L2/3 Exc. <- L2/3 Exc. (n = ', num2str(length(Pre_L23E)), ')'], 'FontWeight', 'normal');
subplot(2, 2, 4); scatter3(Pre_L23I(:, 1), Pre_L23I(:, 2), Pre_L23I(:, 3),...
    10, 'fill', 'MarkerEdgeColor', Clb(2, :), 'MarkerFaceColor', Clb(2, :));
title(['L2/3 Exc. <- L2/3 Inh. (n = ', num2str(length(Pre_L23I)), ')'], 'FontWeight', 'normal');
for k = 1: 4
    subplot(2, 2, k); xlabel('\Deltax (\mum)'); ylabel('\Deltay (\mum)'); zlabel('\Deltaz (\mum)'); view([120 5]);
end
pause(2);
savefig(['Conn_3D_', num2str(-z_boundary), '.fig']);
close;
clear Clb k


%%%%%%%%%%%%% 2D plot of (Δx, Δy) of L4E, L4I, L2/3E, L2/3I presynaptic neurons.
Clb = [0.5 0 0.5; 1 0 1; 1 0 0; 0 0 1];    % L4E, L4I, L2/3E, L2/3I
figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 0.9]);
subplot(2, 2, 1); scatter(Pre_L4E(:, 1), Pre_L4E(:, 2),...
   10, 'fill', 'MarkerEdgeColor', Clb(1, :), 'MarkerFaceColor', Clb(1, :));
title(['L2/3 Exc. <- L4 Exc. (n = ', num2str(length(Pre_L4E)), ')'], 'FontWeight', 'normal');
subplot(2, 2, 2); scatter(Pre_L4I(:, 1), Pre_L4I(:, 2),...
   10, 'fill', 'MarkerEdgeColor', Clb(2, :), 'MarkerFaceColor', Clb(2, :));
title(['L2/3 Exc. <- L4 Inh. (n = ', num2str(length(Pre_L4I)), ')'], 'FontWeight', 'normal');
subplot(2, 2, 3); scatter(Pre_L23E(:, 1), Pre_L23E(:, 2),...
   10, 'fill', 'MarkerEdgeColor', Clb(3, :), 'MarkerFaceColor', Clb(3, :));
title(['L2/3 Exc. <- L2/3 Exc. (n = ', num2str(length(Pre_L23E)), ')'], 'FontWeight', 'normal');
subplot(2, 2, 4); scatter(Pre_L23I(:, 1), Pre_L23I(:, 2),...
   10, 'fill', 'MarkerEdgeColor', Clb(4, :), 'MarkerFaceColor', Clb(4, :));
title(['L2/3 Exc. <- L2/3 Inh. (n = ', num2str(length(Pre_L23I)), ')'], 'FontWeight', 'normal');
for k = 1: 4
   subplot(2, 2, k); axis square; grid on; axis([-400 400 -400 400]);
   set(gca, 'XTick', -400: 100: 400, 'YTick', -400: 100: 400);
   xlabel('\Deltax (\mum)'); ylabel('\Deltay (\mum)');
end
pause(2);
print(gcf, '-dpng', ['Conn_xy_', num2str(-z_boundary), '.png']);
close;
clear Clb k
