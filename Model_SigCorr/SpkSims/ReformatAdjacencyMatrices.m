function [Wrf, Wrr, Kx, KxMax, KxCumsum, Kr, KeMax, KiMax] =...
    ReformatAdjacencyMatrices(dir_JN)

load(dir_JN, 'JN_ef', 'JN_if');
Nf = size(JN_ef, 2); Ne = size(JN_ef, 1); Ni = size(JN_if, 1);
%
Wrf = zeros(sum(JN_ef(:)) + sum(JN_if(:)), 1, 'int32');
Kx = zeros(Nf, 1);
idx_start = 1;
for j = 1: Nf
    Wrf_j = [hist2list(JN_ef(:, j)); hist2list(JN_if(:, j)) + Ne];
    Kx_j = length(Wrf_j); Kx(j) = Kx_j;
    Wrf(idx_start: Kx_j + idx_start - 1) = Wrf_j;
    idx_start = Kx_j + idx_start;
end
clear Wrf_j Kx_j idx_start
%
KxMax = max(Kx);
tmp = cumsum(Kx); KxCumsum = [0; tmp(1: end - 1)]; clear tmp
clear JN_ef JN_if


load(dir_JN, 'JN_ee', 'JN_ie');
Wre = zeros(sum(JN_ee(:)) + sum(JN_ie(:)), 1, 'int32');
Ke = zeros(Ne, 1);
idx_start = 1;
for j = 1: Ne
    Wre_j = [hist2list(JN_ee(:, j)); hist2list(JN_ie(:, j)) + Ne];
    Ke_j = length(Wre_j); Ke(j) = Ke_j;
    Wre(idx_start: Ke_j + idx_start - 1) = Wre_j;
    idx_start = Ke_j + idx_start;
end
clear Wre_j Ke_j idx_start
%
KeMax = max(Ke);
clear JN_ee JN_ie


load(dir_JN, 'JN_ei', 'JN_ii');
Wri = zeros(sum(JN_ei(:)) + sum(JN_ii(:)), 1, 'int32');
Ki = zeros(Ni, 1);
idx_start = 1;
for j = 1: Ni
    Wri_j = [hist2list(JN_ei(:, j)); hist2list(JN_ii(:, j)) + Ne];
    Ki_j = length(Wri_j); Ki(j) = Ki_j;
    Wri(idx_start: Ki_j + idx_start - 1) = Wri_j;
    idx_start = Ki_j + idx_start;
end
clear Wri_j Ki_j idx_start
%
KiMax = max(Ki);
clear JN_ei JN_ii

Kr = [Ke; Ki];

Wrr = [Wre; Wri];



%Kx: (1, Nx), out-degree of each neuon (regardless of post e or i)
%KxCumsum: [0 cumsum(1: end - 1)]
%x x x x v v d d d
%4 2 3
%0 4 6 9
%KxMax

%Ke
%KeMax

%Ki
%KiMax


