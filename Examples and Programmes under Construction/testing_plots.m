%% Comparing the effect of different estimation size (EnKF, DExKF_IFAC:1D,2D)
clear all
close all
clc

disp('--------------------------')
disp('EnKF')
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:2D_IFAC')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:3D_IFAC')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE3(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:4D_IFAC')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_4DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE4(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)

figure, plot(RMSE0), hold on,
plot(RMSE1),
plot(RMSE2),
plot(RMSE3),
plot(RMSE4),
legend('EnKF','DExKF:1D','DExKF:2D','DExKF:3D','DExKF:4D')

%% Comparing different fusing algorithms (EnKF, DExKF_1D:IFAC,CIN,CI,EI,ICI)
clear all
close all
clc

disp('--------------------------')
disp('EnKF')
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_CIN')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_CIN_1DZ\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% disp('--------------------------')
% disp('DExKF:1D_CI')
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_CI_1DZ\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% RMSE        = 0;
% maxError    = 0;
% time        = 0;
% for i = 1:Wp.sim.NN
%     RMSE        = RMSE + sol_array(i).score.RMSE_cline;
%     maxError    = maxError + sol_array(i).score.maxError;
%     time        = time + sol_array(i).score.CPUtime;
%     RMSE3(i) = sol_array(i).score.RMSE_cline;
% end
% avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
% avg_maxError = maxError/(Wp.sim.NN - 1)
% avg_time = time/(Wp.sim.NN - 1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% disp('--------------------------')
% disp('DExKF:1D_EI')
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_EI_1DZ\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% RMSE        = 0;
% maxError    = 0;
% time        = 0;
% for i = 1:Wp.sim.NN
%     RMSE        = RMSE + sol_array(i).score.RMSE_cline;
%     maxError    = maxError + sol_array(i).score.maxError;
%     time        = time + sol_array(i).score.CPUtime;
%     RMSE4(i) = sol_array(i).score.RMSE_cline;
% end
% avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
% avg_maxError = maxError/(Wp.sim.NN - 1)
% avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_ICI')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_ICI_1DZ\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim5,timeWFSim5,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE5(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)

figure, plot(RMSE0), hold on,
plot(RMSE1),
plot(RMSE2),
% plot(RMSE3),
% plot(RMSE4),
plot(RMSE5),
legend('EnKF','DExKF:IFAC','DExKF:CIN','DExKF:ICI')
% legend('EnKF','DExKF:IFAC','DExKF:CIN','DExKF:CI','DExKF:EI','DExKF:ICI')

%% Comparing the effects of linearizing the non-linear system at different frequencies (DExKF:1D_IFAC: 10,20,50,100,Inf)
clear all
close all
clc

disp('--------------------------')
disp('EnKF')
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC: 10')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_10\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC: 20')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_20\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC: 50')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_50\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE3(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC: 100')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_100\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE4(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC: Inf')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim5,timeWFSim5,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE5(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)
avg_maxError = maxError/(Wp.sim.NN - 1)
avg_time = time/(Wp.sim.NN - 1)

figure, plot(RMSE0), hold on,
plot(RMSE1),
plot(RMSE2),
plot(RMSE3),
plot(RMSE4),
plot(RMSE5),
legend('EnKF','DExKF:10','DExKF:20','DExKF:50','DExKF:100','DExKF:Inf')

%% Disturbances added to the initial flow field (EnKF, DExKF:1D_IFAC)
clear all
close all
clc

load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_noise_init_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
for i = 1:Wp.sim.NN
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_noise_init_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    RMSE2(i) = sol_array(i).score.RMSE_cline;
end
figure, plot(RMSE1), hold on, plot(RMSE2)
legend('EnKF','DExKF')

%% Disturbances added to the output (EnKF, DExKF:1D_IFAC)
clear all
close all
clc

load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_noise_obs_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
for i = 1:Wp.sim.NN
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
% plot(u), hold on, legend('EnKF', 'DExKF')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_noise_obs_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    RMSE2(i) = sol_array(i).score.RMSE_cline;
end
figure, plot(RMSE1), hold on, plot(RMSE2)
legend('EnKF','DExKF')
%% Estimate the freestream flow (EnKF, DExKF:1D_IFAC)
clear all
close all
clc

load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_uinf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
u1 = [];
for i = 1:Wp.sim.NN
    u1(i) = sol_array(i).site.u_Inf;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
% plot(u), hold on, legend('EnKF', 'DExKF')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_uinf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
u2 = [];
for i = 1:Wp.sim.NN
    u2(i) = sol_array(i).site.u_Inf;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
end
figure, plot(u1), hold on, 
plot(u2), plot(Wp.site.actual_u_Inf*ones(1,Wp.sim.NN)) 
legend('EnKF','DExKF'), axis([0 2000 0 9])
figure, plot(RMSE1), hold on, plot(RMSE2)
legend('EnKF','DExKF')
%%
clear all
close all
clc
% % load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/axi_2turb_alm_turb_ExKF/ExKF_2turb_workspace.mat')
% load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_ExKF\ExKF_2turb_workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim1,timeWFSim1,pwLES,timeLES] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/axi_2turb_alm_turb_EnKF/EnKF_2turb_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_EnKF\EnKF_2turb_workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% % load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 2/dexkf_2turb_IFAC_1D_workspace.mat')
% load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 2\dexkf_2turb_IFAC_1D_workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 3/dexkf_2turb_IFAC_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 3\dexkf_2turb_IFAC_1D_workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 4/dexkf_2turb_IFAC_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 4\dexkf_2turb_IFAC_1D_workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim5,timeWFSim5,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% % load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC2DZ/trail 2/dexkf_2turb_IFAC_2D_workspace.mat')
% load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC2DZ\trail 2\dexkf_2turb_IFAC_2D_workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim6,timeWFSim6,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_dexkf_ICI1DZ\trail 2\dexkf_2turb_ICI_1D_workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim7,timeWFSim7,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
hold off

% e11 = pwWFSim1(1,:) - pwLES(1,:);
% e12 = pwWFSim1(2,:) - pwLES(2,:);
% 
% e21 = pwWFSim2(1,:) - pwLES(1,:);
% e22 = pwWFSim2(2,:) - pwLES(2,:);
% 
% e31 = pwWFSim3(1,:) - pwLES(1,:);
% e32 = pwWFSim3(2,:) - pwLES(2,:);
% 
% e41 = pwWFSim4(1,:) - pwLES(1,:);
% e42 = pwWFSim4(2,:) - pwLES(2,:);
% 
% e51 = pwWFSim5(1,:) - pwLES(1,:);
% e52 = pwWFSim5(2,:) - pwLES(2,:);
% 
% figure,plot(timeWFSim1,e11), hold on, plot(timeWFSim1,e21),plot(timeWFSim1,e31), plot(timeWFSim1,e41), hold off
% figure,plot(timeWFSim1,e12), hold on, plot(timeWFSim1,e22),plot(timeWFSim1,e32), plot(timeWFSim1,e42), hold off
% 
% ee11 = sum(e11)/length(e11)
% ee21 = sum(e21)/length(e21)
% ee31 = sum(e31)/length(e31)
% ee41 = sum(e41)/length(e41)
% 
% ee12 = sum(e12)/length(e12)
% ee22 = sum(e22)/length(e22)
% ee32 = sum(e32)/length(e32)
% ee42 = sum(e42)/length(e42)

%% 2Turbine

clear all
close all
clc

disp('--------------------------')
disp('ExKF: Linearize every iterations')
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/axi_2turb_alm_turb_ExKF/ExKF_2turb_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_ExKF\ExKF_2turb_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('EnKF')
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/axi_2turb_alm_turb_EnKF/EnKF_2turb_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_EnKF\EnKF_2turb_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF: 1D + IFAC: Trail 1 (Unoptimized + linearize once every 5 iterations)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 1/dexkf_2turb_IFAC_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 1\dexkf_2turb_IFAC_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF: 1D + IFAC: Trail 2 (Optimized + linearize once every 50 iterations)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 2/dexkf_2turb_IFAC_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 2\dexkf_2turb_IFAC_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF: 1D + IFAC: Trail 3 (Optimized + linearize once every 5 iterations)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 3/dexkf_2turb_IFAC_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 3\dexkf_2turb_IFAC_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF: 1D + IFAC: Trail 4 (Optimized + linearize at k = 1 + Conv. KF)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC1DZ/trail 3/dexkf_2turb_IFAC_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC1DZ\trail 4\dexkf_2turb_IFAC_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('DExKF: 2D + IFAC')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC2DZ/trail 1/dexkf_2turb_IFAC_2D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC2DZ\trail 1\dexkf_2turb_IFAC_2D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF: 2D + IFAC: Trail 2 (Optimized + linearize once every 50 iterations)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC2DZ/trail 2/dexkf_2turb_IFAC_2D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC2DZ\trail 2\dexkf_2turb_IFAC_2D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('DExKF: 3D + IFAC')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC3DZ/dexkf_2turb_IFAC_3D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC3DZ\dexkf_2turb_IFAC_3D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('DExKF: 4D + IFAC')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/IFAC/axi_2turb_alm_turb_dexkf_IFAC4DZ/dexkf_2turb_IFAC_4D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\IFAC\axi_2turb_alm_turb_dexkf_IFAC4DZ\dexkf_2turb_IFAC_4D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('DExKF: 1D + ICI: Trail 1 (Unoptimized + linearize once every 5 iterations)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/axi_2turb_alm_turb_dexkf_ICI1DZ/dexkf_2turb_ICI_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_dexkf_ICI1DZ\trail 1\dexkf_2turb_ICI_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('DExKF: 1D + ICI: Trail 2 (Optimized + linearize at k = 1)')
clear all
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/2Turbine/axi_2turb_alm_turb_dexkf_ICI1DZ/dexkf_2turb_ICI_1D_workspace.mat')
load('H:\MSc Thesis\WFOBS - Queue\2Turbine\axi_2turb_alm_turb_dexkf_ICI1DZ\trail 2\dexkf_2turb_ICI_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%% 9Turbine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('9TURBINE: SIM')
clear all
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/9Turbine/apc_9turb_alm_turb/apc_9turb_alm_turb_sim/9turb_sim_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------')
disp('9TURBINE: DExKF: 1D + IFAC')
clear all
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/9Turbine/apc_9turb_alm_turb/apc_9turb_alm_turb_dexkf_IFAC1DZ/Trail 1/dexkf_9turb_IFAC_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

disp('--------------------------')
disp('9TURBINE: DExKF: 1D + IFAC')
clear all
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/9Turbine/apc_9turb_alm_turb/apc_9turb_alm_turb_dexkf_IFAC1DZ/Trail 2/dexkf_9turb_IFAC_1D_workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );

disp('--------------------------')
disp('9TURBINE (apc_9turb_adm_noturb): DExKF: 1D + IFAC')
clear all
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFOBS - Queue/9Turbine/apc_9turb_adm_noturb/apc_9turb_adm_noturb_dexkf_IFAC1DZ/workspace.mat')
RMSE = 0;
% figure, hold on,
for i = 1:Wp.sim.NN
RMSE = RMSE + sol_array(i).score.RMSE_cline;
% plot(i,sol_array(i).score.RMSE_cline,'o')
end
% hold off
% RMSE
avg_RMSE_cline = RMSE/(Wp.sim.NN - 1)

maxError = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
maxError = maxError + sol_array(i).score.maxError;
% plot(i,sol_array(i).score.maxError,'o')
end
% hold off
% maxError
avg_maxError = maxError/(Wp.sim.NN - 1)

time = 0;
% figure, hold on,
for i = 2:Wp.sim.NN
time = time + sol_array(i).score.CPUtime;
% plot(i,sol_array(i).score.CPUtime,'o')
end
% hold off
% time
avg_time = time/(Wp.sim.NN - 1)
% plotWFObs( Wp,sol_array,sys,scriptOptions,strucObs );