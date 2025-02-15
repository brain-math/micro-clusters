function Parameters_Recurrent_dthetadep(Parset_i_spatial, Parset_i_dthetadep, ifDthetaDep, sigma_dtheta, a_dtheta)

dir_JN_spatial = [pwd, '/Parameters/JN_', num2str(Parset_i_spatial), '.mat'];
dir_JN_dthetadep = [pwd, '/Parameters/JN_', num2str(Parset_i_dthetadep), '.mat'];
%
load([pwd, '/Results/Results_Total_', num2str(Parset_i_spatial), '.mat'],...
    'd_tot', 'Pref_theta_F', 'Pref_theta_E', 'Pref_theta_I');
load([pwd, '/Parameters/Parameters_Recurrent_', num2str(Parset_i_spatial), '.mat'], 'param');
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


