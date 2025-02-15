function Parameters_Recurrent_dthetadep(Parset_i_spatial, Parset_i_dthetadep, ifDthetaDep, sigma_dtheta, a_dtheta)

rng('shuffle');


if Parset_i_spatial == 4
    dir_JN_spatial = [pwd, '/Parameters/JN_1.mat'];
else
    dir_JN_spatial = [pwd, '/Parameters/JN_', num2str(Parset_i_spatial), '.mat'];
end
dir_JN_dthetadep = [pwd, '/Parameters/JN_', num2str(Parset_i_dthetadep), '.mat'];
%
load([pwd, '/Parameters/Parameters_FFWD.mat'], 'd_tot');
load([pwd, '/Results/Results_Total_', num2str(Parset_i_spatial), '.mat'],...
    'Pref_theta_F', 'Pref_theta_E', 'Pref_theta_I');
if Parset_i_spatial == 4
    load([pwd, '/Parameters/Parameters_Recurrent_1.mat'], 'param');
else
    load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i_spatial), '.mat'], 'param');
end
rmfield(param, 'Kx'); rmfield(param, 'KxMax'); rmfield(param, 'KxCumsum'); 
rmfield(param, 'Kr'); rmfield(param, 'KeMax'); rmfield(param, 'KiMax'); 

param.ifDthetaDep = ifDthetaDep;
param.sigma_dtheta = sigma_dtheta;
param.a_dtheta = a_dtheta;

getRecDthetaDepAdjacencyMatrices(dir_JN_spatial,...
    Pref_theta_F, Pref_theta_E, Pref_theta_I, d_tot,...
    ifDthetaDep, sigma_dtheta, a_dtheta, dir_JN_dthetadep);


[Wrf, Wrr, param.Kx, param.KxMax, param.KxCumsum,...
    param.Kr, param.KeMax, param.KiMax] =...
    ReformatAdjacencyMatrices(dir_JN_dthetadep);


save([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i_dthetadep), '.mat'], 'param', 'Wrf', 'Wrr');
fprintf(['Parameters_Recurrent_', num2str(Parset_i_dthetadep), '.mat saved.\n']);


