%% Different initial value for mixing length and it is not estimated (EnKF, DExKF_IFAC:1D,2D)
% With E
clear all
close all
clc
%SIM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_sim_lmu_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
for i = 1:Wp.sim.NN
    u0(i) = sol_array(i).site.u_Inf;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
    maxError0(i) = sol_array(i).score.maxError;
    RMSE_flow0(i) = sol_array(i).score.RMSE_flow;
end
Wp0 = Wp; sol_array0 = sol_array; sys0 = sys; 
scriptOptions0 = scriptOptions; strucObs0 = strucObs;
%EnKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_lmu_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u1(i) = sol_array(i).site.u_Inf;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
    maxError1(i) = sol_array(i).score.maxError;
    RMSE_flow1(i) = sol_array(i).score.RMSE_flow;
end
Wp1 = Wp; sol_array1 = sol_array; sys1 = sys; 
scriptOptions1 = scriptOptions; strucObs1 = strucObs;
%IFAC_1DZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZE_lmu_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u2(i) = sol_array(i).site.u_Inf;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
    maxError2(i) = sol_array(i).score.maxError;
    RMSE_flow2(i) = sol_array(i).score.RMSE_flow;
end
Wp2 = Wp; sol_array2 = sol_array; sys2 = sys; 
scriptOptions2 = scriptOptions; strucObs2 = strucObs;
%IFAC_2DZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_lmu_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim11,timeWFSim11,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u3(i) = sol_array(i).site.u_Inf;
    RMSE3(i) = sol_array(i).score.RMSE_cline;
    maxError3(i) = sol_array(i).score.maxError;
    RMSE_flow3(i) = sol_array(i).score.RMSE_flow;
end
Wp3 = Wp; sol_array3 = sol_array; sys3 = sys; 
scriptOptions3 = scriptOptions; strucObs3 = strucObs;
%IFAC_3DZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZE_lmu_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim12,timeWFSim12,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u4(i) = sol_array(i).site.u_Inf;
    RMSE4(i) = sol_array(i).score.RMSE_cline;
    maxError4(i) = sol_array(i).score.maxError;
    RMSE_flow4(i) = sol_array(i).score.RMSE_flow;
end
Wp4 = Wp; sol_array4 = sol_array; sys4 = sys; 
scriptOptions4 = scriptOptions; strucObs4 = strucObs;

figure, plot(RMSE0), hold on, 
plot(RMSE1),plot(RMSE2), plot(RMSE3), plot(RMSE4) 
legend('OPEN','EnKF','DExKF_IFAC:1D','DExKF_IFAC:2D','DExKF_IFAC:3D')
title('RMSE_cline'), xlabel('Time (sec)'), ylabel('RMSE')
figure, plot(maxError0), hold on, 
plot(maxError1),plot(maxError2), plot(maxError3), plot(maxError4)
legend('OPEN','EnKF','DExKF_IFAC:1D','DExKF_IFAC:2D','DExKF_IFAC:3D')
title('RMSE'), xlabel('Time (sec)'), ylabel('maxError')
figure, plot(RMSE_flow0), hold on, 
plot(RMSE_flow1),plot(RMSE_flow2), plot(RMSE_flow3), plot(RMSE_flow4)
legend('OPEN','EnKF','DExKF_IFAC:1D','DExKF_IFAC:2D','DExKF_IFAC:3D')
title('RMSE'), xlabel('Time (sec)'), ylabel('RMSE_flow')

hold off
plotWFObs( Wp0,sol_array0,sys0,scriptOptions0,strucObs0 );
plotWFObs( Wp1,sol_array1,sys1,scriptOptions1,strucObs1 );
plotWFObs( Wp2,sol_array2,sys2,scriptOptions2,strucObs2 );
plotWFObs( Wp3,sol_array3,sys3,scriptOptions3,strucObs3 );
plotWFObs( Wp4,sol_array4,sys4,scriptOptions4,strucObs4 );
%% Noise
clear all
close all
clc
disp('--------------------------')
disp('DExKF:1D_IFAC')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_NL_Inf\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
RMSE        = 0;
maxError    = 0;
time        = 0;
for i = 1:Wp.sim.NN
    RMSE        = RMSE + sol_array(i).score.RMSE_cline;
    maxError    = maxError + sol_array(i).score.maxError;
    time        = time + sol_array(i).score.CPUtime;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
end
avg_RMSE_cline1 = RMSE/(Wp.sim.NN - 1)
avg_maxError1 = maxError/(Wp.sim.NN - 1)
avg_time1 = time/(Wp.sim.NN - 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
disp('--------------------------')
disp('DExKF:1D_IFAC')
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_noise1\workspace.mat')
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
avg_RMSE_cline2 = RMSE/(Wp.sim.NN - 1)
avg_maxError2 = maxError/(Wp.sim.NN - 1)
avg_time2 = time/(Wp.sim.NN - 1)

figure, plot(RMSE1), hold on,
plot(RMSE2),
legend('DExKF_1D: No Noise','DExKF_1D: Noise')

%% Different initial value for freestream velocity and it is not estimated (EnKF, DExKF_IFAC:1D,2D)
% Without E
clear all
close all
clc
%SIM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_sim_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
for i = 1:Wp.sim.NN
    u0(i) = sol_array(i).site.u_Inf;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
    maxError0(i) = sol_array(i).score.maxError;
    RMSE_flow0(i) = sol_array(i).score.RMSE_flow;
end
Wp0 = Wp; sol_array0 = sol_array; sys0 = sys; 
scriptOptions0 = scriptOptions; strucObs0 = strucObs;
%EnKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u1(i) = sol_array(i).site.u_Inf;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
    maxError1(i) = sol_array(i).score.maxError;
    RMSE_flow1(i) = sol_array(i).score.RMSE_flow;
end
Wp1 = Wp; sol_array1 = sol_array; sys1 = sys; 
scriptOptions1 = scriptOptions; strucObs1 = strucObs;
%ExKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_ExKF_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim2,timeWFSim2,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u2(i) = sol_array(i).site.u_Inf;
    RMSE2(i) = sol_array(i).score.RMSE_cline;
    maxError2(i) = sol_array(i).score.maxError;
    RMSE_flow2(i) = sol_array(i).score.RMSE_flow;
end
Wp2 = Wp; sol_array2 = sol_array; sys2 = sys; 
scriptOptions2 = scriptOptions; strucObs2 = strucObs;
%IFAC_1DZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u3(i) = sol_array(i).site.u_Inf;
    RMSE3(i) = sol_array(i).score.RMSE_cline;
    maxError3(i) = sol_array(i).score.maxError;
    RMSE_flow3(i) = sol_array(i).score.RMSE_flow;
end
Wp3 = Wp; sol_array3 = sol_array; sys3 = sys; 
scriptOptions3 = scriptOptions; strucObs3 = strucObs;
%IFAC_2DZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZ_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u4(i) = sol_array(i).site.u_Inf;
    RMSE4(i) = sol_array(i).score.RMSE_cline;
    maxError4(i) = sol_array(i).score.maxError;
    RMSE_flow4(i) = sol_array(i).score.RMSE_flow;
end
Wp4 = Wp; sol_array4 = sol_array; sys4 = sys; 
scriptOptions4 = scriptOptions; strucObs4 = strucObs;
%IFAC_4DZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_4DZ_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim5,timeWFSim5,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u5(i) = sol_array(i).site.u_Inf;
    RMSE5(i) = sol_array(i).score.RMSE_cline;
    maxError5(i) = sol_array(i).score.maxError;
    RMSE_flow5(i) = sol_array(i).score.RMSE_flow;
end
Wp5 = Wp; sol_array5 = sol_array; sys5 = sys; 
scriptOptions5 = scriptOptions; strucObs5 = strucObs;
%CIN_4DZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_CIN_4DZ_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim6,timeWFSim6,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u6(i) = sol_array(i).site.u_Inf;
    RMSE6(i) = sol_array(i).score.RMSE_cline;
    maxError6(i) = sol_array(i).score.maxError;
    RMSE_flow6(i) = sol_array(i).score.RMSE_flow;
end
Wp6 = Wp; sol_array6 = sol_array; sys6 = sys; 
scriptOptions6 = scriptOptions; strucObs6 = strucObs;
%CIN_1DC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_CIN_1DC_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim7,timeWFSim7,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u7(i) = sol_array(i).site.u_Inf;
    RMSE7(i) = sol_array(i).score.RMSE_cline;
    maxError7(i) = sol_array(i).score.maxError;
    RMSE_flow7(i) = sol_array(i).score.RMSE_flow;
end
Wp7 = Wp; sol_array7 = sol_array; sys7 = sys; 
scriptOptions7 = scriptOptions; strucObs7 = strucObs;
%CIN_2DC%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_CIN_2DC_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim8,timeWFSim8,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u8(i) = sol_array(i).site.u_Inf;
    RMSE8(i) = sol_array(i).score.RMSE_cline;
    maxError8(i) = sol_array(i).score.maxError;
    RMSE_flow8(i) = sol_array(i).score.RMSE_flow;
end
Wp8 = Wp; sol_array8 = sol_array; sys8 = sys; 
scriptOptions8 = scriptOptions; strucObs8 = strucObs;
%CIN_1DZ%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_CIN_1DZ_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim9,timeWFSim9,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u9(i) = sol_array(i).site.u_Inf;
    RMSE9(i) = sol_array(i).score.RMSE_cline;
    maxError9(i) = sol_array(i).score.maxError;
    RMSE_flow9(i) = sol_array(i).score.RMSE_flow;
end
Wp9 = Wp; sol_array9 = sol_array; sys9 = sys; 
scriptOptions9 = scriptOptions; strucObs9 = strucObs;
figure, plot(RMSE0), hold on, 
plot(RMSE1),plot(RMSE2), plot(RMSE3), plot(RMSE4), plot(RMSE5), 
plot(RMSE6), plot(RMSE7), plot(RMSE8), plot(RMSE9), 
legend('OPEN','EnKF','ExKF','DExKF_IFAC:1D','DExKF_IFAC:2D','DExKF_IFAC:4D','DExKF_CIN:4DZ','DExKF_CIN:1DC','DExKF_CIN:2DC','DExKF_CIN:1DZ')
title('RMSE_cline'), xlabel('Time (sec)'), ylabel('RMSE')
figure, plot(maxError0), hold on, 
plot(maxError1),plot(maxError2), plot(maxError3), plot(maxError4), plot(maxError5), 
plot(maxError6), plot(maxError7), plot(maxError8), plot(maxError9), 
legend('OPEN','EnKF','ExKF','DExKF_IFAC:1D','DExKF_IFAC:2D','DExKF_IFAC:4D','DExKF_CIN:4DZ','DExKF_CIN:1DC','DExKF_CIN:2DC','DExKF_CIN:1DZ')
title('RMSE'), xlabel('Time (sec)'), ylabel('maxError')
figure, plot(RMSE_flow0), hold on, 
plot(RMSE_flow1),plot(RMSE_flow2), plot(RMSE_flow3), plot(RMSE_flow4), plot(RMSE_flow5), 
plot(RMSE_flow6), plot(RMSE_flow7), plot(RMSE_flow8), plot(RMSE_flow9), 
legend('OPEN','EnKF','ExKF','DExKF_IFAC:1D','DExKF_IFAC:2D','DExKF_IFAC:4D','DExKF_CIN:4DZ','DExKF_CIN:1DC','DExKF_CIN:2DC','DExKF_CIN:1DZ')
title('RMSE'), xlabel('Time (sec)'), ylabel('RMSE_flow')

hold off
plotWFObs( Wp0,sol_array0,sys0,scriptOptions0,strucObs0 );
plotWFObs( Wp1,sol_array1,sys1,scriptOptions1,strucObs1 );
plotWFObs( Wp2,sol_array2,sys2,scriptOptions2,strucObs2 );
plotWFObs( Wp3,sol_array3,sys3,scriptOptions3,strucObs3 );
plotWFObs( Wp4,sol_array4,sys4,scriptOptions4,strucObs4 );
plotWFObs( Wp5,sol_array5,sys5,scriptOptions5,strucObs5 );
plotWFObs( Wp6,sol_array6,sys6,scriptOptions6,strucObs6 );
plotWFObs( Wp7,sol_array7,sys7,scriptOptions7,strucObs7 );
plotWFObs( Wp8,sol_array8,sys8,scriptOptions8,strucObs8 );
plotWFObs( Wp9,sol_array9,sys9,scriptOptions9,strucObs9 );
%% Different initial value for freestream velocity and it is not estimated (EnKF, DExKF_IFAC:1D,2D)
% With E
clear all
close all
clc
%SIM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_sim_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim0,timeWFSim0,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 1 );
for i = 1:Wp.sim.NN
    u0(i) = sol_array(i).site.u_Inf;
    RMSE0(i) = sol_array(i).score.RMSE_cline;
    maxError0(i) = sol_array(i).score.maxError;
    RMSE_flow0(i) = sol_array(i).score.RMSE_flow;
end
Wp0 = Wp; sol_array0 = sol_array; sys0 = sys; 
scriptOptions0 = scriptOptions; strucObs0 = strucObs;
%EnKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_EnKF_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim1,timeWFSim1,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u1(i) = sol_array(i).site.u_Inf;
    RMSE1(i) = sol_array(i).score.RMSE_cline;
    maxError1(i) = sol_array(i).score.maxError;
    RMSE_flow1(i) = sol_array(i).score.RMSE_flow;
end
Wp1 = Wp; sol_array1 = sol_array; sys1 = sys; 
scriptOptions1 = scriptOptions; strucObs1 = strucObs;
% %ExKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_ExKF_uinf_noest\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     u4(i) = sol_array(i).site.u_Inf;
%     RMSE4(i) = sol_array(i).score.RMSE_cline;
%     maxError4(i) = sol_array(i).score.maxError;
%     RMSE_flow4(i) = sol_array(i).score.RMSE_flow;
% end
% Wp4 = Wp; sol_array4 = sol_array; sys4 = sys; 
% scriptOptions4 = scriptOptions; strucObs4 = strucObs;
% %IFAC_1DZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZE_uinf_noest\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim10,timeWFSim10,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     u10(i) = sol_array(i).site.u_Inf;
%     RMSE10(i) = sol_array(i).score.RMSE_cline;
%     maxError10(i) = sol_array(i).score.maxError;
%     RMSE_flow10(i) = sol_array(i).score.RMSE_flow;
% end
% Wp10 = Wp; sol_array10 = sol_array; sys10 = sys; 
% scriptOptions10 = scriptOptions; strucObs10 = strucObs;
%IFAC_2DZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim11,timeWFSim11,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u11(i) = sol_array(i).site.u_Inf;
    RMSE11(i) = sol_array(i).score.RMSE_cline;
    maxError11(i) = sol_array(i).score.maxError;
    RMSE_flow11(i) = sol_array(i).score.RMSE_flow;
end
Wp11 = Wp; sol_array11 = sol_array; sys11 = sys; 
scriptOptions11 = scriptOptions; strucObs11 = strucObs;
%IFAC_3DZE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZE_uinf_noest\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim12,timeWFSim12,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u12(i) = sol_array(i).site.u_Inf;
    RMSE12(i) = sol_array(i).score.RMSE_cline;
    maxError12(i) = sol_array(i).score.maxError;
    RMSE_flow12(i) = sol_array(i).score.RMSE_flow;
end
Wp12 = Wp; sol_array12 = sol_array; sys12 = sys; 
scriptOptions12 = scriptOptions; strucObs12 = strucObs;
% %NoFusion_1DCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_nofusion_1DC_uinf_noest\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim13,timeWFSim13,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     u13(i) = sol_array(i).site.u_Inf;
%     RMSE13(i) = sol_array(i).score.RMSE_cline;
%     maxError13(i) = sol_array(i).score.maxError;
%     RMSE_flow13(i) = sol_array(i).score.RMSE_flow;
% end
% Wp13 = Wp; sol_array13 = sol_array; sys13 = sys; 
% scriptOptions13 = scriptOptions; strucObs13 = strucObs;
% %NoFusion_2DCE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_nofusion_2DC_uinf_noest\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim14,timeWFSim14,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     u14(i) = sol_array(i).site.u_Inf;
%     RMSE14(i) = sol_array(i).score.RMSE_cline;
%     maxError14(i) = sol_array(i).score.maxError;
%     RMSE_flow14(i) = sol_array(i).score.RMSE_flow;
% end
% Wp14 = Wp; sol_array14 = sol_array; sys14 = sys; 
% scriptOptions14 = scriptOptions; strucObs14 = strucObs;
% %IFAC_1DZE_opti%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZE_uinf_noest_opti\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim15,timeWFSim15,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     u15(i) = sol_array(i).site.u_Inf;
%     RMSE15(i) = sol_array(i).score.RMSE_cline;
%     maxError15(i) = sol_array(i).score.maxError;
%     RMSE_flow15(i) = sol_array(i).score.RMSE_flow;
% end
% Wp15 = Wp; sol_array15 = sol_array; sys15 = sys; 
% scriptOptions15 = scriptOptions; strucObs15 = strucObs;
% %IFAC_2DZE_opti%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_opti\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim16,timeWFSim16,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     u16(i) = sol_array(i).site.u_Inf;
%     RMSE16(i) = sol_array(i).score.RMSE_cline;
%     maxError16(i) = sol_array(i).score.maxError;
%     RMSE_flow16(i) = sol_array(i).score.RMSE_flow;
% end
% Wp16 = Wp; sol_array16 = sol_array; sys16 = sys; 
% scriptOptions16 = scriptOptions; strucObs16 = strucObs;
% %IFAC_3DZE_opti%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% clear scriptOptions sol_array strucObs sys Wp
% load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZE_uinf_noest_opti\workspace.mat')
% scriptOptions.powerForecast = 0;
% [pwWFSim17,timeWFSim17,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
% for i = 1:Wp.sim.NN
%     u17(i) = sol_array(i).site.u_Inf;
%     RMSE17(i) = sol_array(i).score.RMSE_cline;
%     maxError17(i) = sol_array(i).score.maxError;
%     RMSE_flow17(i) = sol_array(i).score.RMSE_flow;
% end
% Wp17 = Wp; sol_array17 = sol_array; sys17 = sys; 
% scriptOptions17 = scriptOptions; strucObs17 = strucObs;
%IFAC_2DZE_superopti_1em2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_superopti_1em2\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim18,timeWFSim18,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u18(i) = sol_array(i).site.u_Inf;
    RMSE18(i) = sol_array(i).score.RMSE_cline;
    maxError18(i) = sol_array(i).score.maxError;
    RMSE_flow18(i) = sol_array(i).score.RMSE_flow;
end
Wp18 = Wp; sol_array18 = sol_array; sys18 = sys; 
scriptOptions18 = scriptOptions; strucObs18 = strucObs;
%IFAC_2DZE_superopti_1em3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_superopti_1em3\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim19,timeWFSim19,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u19(i) = sol_array(i).site.u_Inf;
    RMSE19(i) = sol_array(i).score.RMSE_cline;
    maxError19(i) = sol_array(i).score.maxError;
    RMSE_flow19(i) = sol_array(i).score.RMSE_flow;
end
Wp19 = Wp; sol_array19 = sol_array; sys19 = sys; 
scriptOptions19 = scriptOptions; strucObs19 = strucObs;
%IFAC_2DZE_superopti_1em4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZE_uinf_noest_superopti_1em4\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim20,timeWFSim20,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u20(i) = sol_array(i).site.u_Inf;
    RMSE20(i) = sol_array(i).score.RMSE_cline;
    maxError20(i) = sol_array(i).score.maxError;
    RMSE_flow20(i) = sol_array(i).score.RMSE_flow;
end
Wp20 = Wp; sol_array20 = sol_array; sys20 = sys; 
scriptOptions20 = scriptOptions; strucObs20 = strucObs;
%IFAC_3DZE_superopti_1em4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZE_uinf_noest_superopti_1em4\workspace.mat')
scriptOptions.powerForecast = 0;
% [pwWFSim21,timeWFSim21,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    u21(i) = sol_array(i).site.u_Inf;
    RMSE21(i) = sol_array(i).score.RMSE_cline;
    maxError21(i) = sol_array(i).score.maxError;
    RMSE_flow21(i) = sol_array(i).score.RMSE_flow;
end
Wp21 = Wp; sol_array21 = sol_array; sys21 = sys; 
scriptOptions21 = scriptOptions; strucObs21 = strucObs;

figure, plot(RMSE0), hold on, plot(RMSE1), %plot(RMSE4), plot(RMSE10),
plot(RMSE11),plot(RMSE12),
% plot(RMSE15),plot(RMSE16),plot(RMSE17)
plot(RMSE18),plot(RMSE19),plot(RMSE20),plot(RMSE21)
legend('OPEN','EnKF','DExKF_IFAC:2DE','DExKF_IFAC:3DE','DExKF_IFAC:2DE_superopti_1em2','DExKF_IFAC:2DE_superopti_1em3','DExKF_IFAC:2DE_superopti_1em4','DExKF_IFAC:3DE_superopti_1em4')
title('RMSE_cline'), xlabel('Time (sec)'), ylabel('RMSE')
figure, plot(maxError0), hold on, plot(maxError1), %plot(maxError4), plot(maxError10),
plot(maxError11),plot(maxError12), 
% plot(maxError15),plot(maxError16),plot(maxError17)
plot(maxError18),plot(maxError19),plot(maxError20),plot(maxError21)
legend('OPEN','EnKF','DExKF_IFAC:2DE','DExKF_IFAC:3DE','DExKF_IFAC:2DE_superopti_1em2','DExKF_IFAC:2DE_superopti_1em3','DExKF_IFAC:2DE_superopti_1em4','DExKF_IFAC:3DE_superopti_1em4')
title('maxError'), xlabel('Time (sec)'), ylabel('maxError')
figure, plot(RMSE_flow0), hold on, plot(RMSE_flow1), %plot(RMSE_flow4),plot(RMSE_flow10),
plot(RMSE_flow11),plot(RMSE_flow12), 
% plot(RMSE_flow15),plot(RMSE_flow16),plot(RMSE_flow17)
plot(RMSE_flow18),plot(RMSE_flow19),plot(RMSE_flow20),plot(RMSE_flow21)
legend('OPEN','EnKF','DExKF_IFAC:2DE','DExKF_IFAC:3DE','DExKF_IFAC:2DE_superopti_1em2','DExKF_IFAC:2DE_superopti_1em3','DExKF_IFAC:2DE_superopti_1em4','DExKF_IFAC:3DE_superopti_1em4')
title('RMSE_flow'), xlabel('Time (sec)'), ylabel('RMSE_flow')

hold off
plotWFObs( Wp0,sol_array0,sys0,scriptOptions0,strucObs0 );
plotWFObs( Wp1,sol_array1,sys1,scriptOptions1,strucObs1 );
% % plotWFObs( Wp4,sol_array4,sys4,scriptOptions4,strucObs4 );
% % plotWFObs( Wp10,sol_array10,sys10,scriptOptions10,strucObs10 );
plotWFObs( Wp11,sol_array11,sys11,scriptOptions11,strucObs11 );
plotWFObs( Wp12,sol_array12,sys12,scriptOptions12,strucObs12 );
% % plotWFObs( Wp15,sol_array15,sys15,scriptOptions15,strucObs15 );
% % plotWFObs( Wp16,sol_array16,sys16,scriptOptions16,strucObs16 );
% % plotWFObs( Wp17,sol_array17,sys17,scriptOptions17,strucObs17 );
plotWFObs( Wp18,sol_array18,sys18,scriptOptions18,strucObs18 );
plotWFObs( Wp19,sol_array19,sys19,scriptOptions19,strucObs19 );
plotWFObs( Wp20,sol_array20,sys20,scriptOptions20,strucObs20 );
plotWFObs( Wp21,sol_array21,sys21,scriptOptions21,strucObs21 );
%% Comparing the effect of different estimation size (EnKF, DExKF_IFAC:1D,2D,3D,4D)
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
avg_RMSE_cline0 = RMSE/(Wp.sim.NN - 1)
avg_maxError0 = maxError/(Wp.sim.NN - 1)
avg_time0 = time/(Wp.sim.NN - 1)
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
avg_RMSE_cline1 = RMSE/(Wp.sim.NN - 1)
avg_maxError1 = maxError/(Wp.sim.NN - 1)
avg_time1 = time/(Wp.sim.NN - 1)
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
avg_RMSE_cline2 = RMSE/(Wp.sim.NN - 1)
avg_maxError2 = maxError/(Wp.sim.NN - 1)
avg_time2 = time/(Wp.sim.NN - 1)
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
avg_RMSE_cline3 = RMSE/(Wp.sim.NN - 1)
avg_maxError3 = maxError/(Wp.sim.NN - 1)
avg_time3 = time/(Wp.sim.NN - 1)
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
avg_RMSE_cline4 = RMSE/(Wp.sim.NN - 1)
avg_maxError4 = maxError/(Wp.sim.NN - 1)
avg_time4 = time/(Wp.sim.NN - 1)

figure, plot(RMSE0), hold on,
plot(RMSE1),
plot(RMSE2),
plot(RMSE3),
plot(RMSE4),
legend('EnKF','DExKF:1D','DExKF:2D','DExKF:3D','DExKF:4D')
D = Wp.turbine.Drotor;
dia1    = D;
dia2    = 2*D;
dia3    = 3*D;
dia4    = 4*D;
dia     = [dia1,dia2,dia3,dia4];
avg_time= [avg_time1,avg_time2,avg_time3,avg_time4];
avg_RMSE_cline = [avg_RMSE_cline1,avg_RMSE_cline2,avg_RMSE_cline3,avg_RMSE_cline4];
figure, 
yyaxis left, plot(dia,avg_time,'-.*')
ylabel('Average time (sec)')
yyaxis right,plot(dia,avg_RMSE_cline,':o')
axis([D-(D/10) 4*D+(D/10), 0.3 0.5])
ylabel('Average RMSE Cline (sec)');
xlabel('Radius (m)')
xticks(dia)
xticklabels({'1D','2D','3D','4D'})
title('Effect of different estimation size')

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
avg_RMSE_cline0 = RMSE/(Wp.sim.NN - 1)
avg_maxError0 = maxError/(Wp.sim.NN - 1)
avg_time0 = time/(Wp.sim.NN - 1)
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
avg_RMSE_cline1 = RMSE/(Wp.sim.NN - 1)
avg_maxError1 = maxError/(Wp.sim.NN - 1)
avg_time1 = time/(Wp.sim.NN - 1)
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
avg_RMSE_cline2 = RMSE/(Wp.sim.NN - 1)
avg_maxError2 = maxError/(Wp.sim.NN - 1)
avg_time2 = time/(Wp.sim.NN - 1)
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
% avg_RMSE_cline3 = RMSE/(Wp.sim.NN - 1)
% avg_maxError3 = maxError/(Wp.sim.NN - 1)
% avg_time3 = time/(Wp.sim.NN - 1)
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
% avg_RMSE_cline4 = RMSE/(Wp.sim.NN - 1)
% avg_maxError4 = maxError/(Wp.sim.NN - 1)
% avg_time4 = time/(Wp.sim.NN - 1)
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
avg_RMSE_cline5 = RMSE/(Wp.sim.NN - 1)
avg_maxError5 = maxError/(Wp.sim.NN - 1)
avg_time5 = time/(Wp.sim.NN - 1)

figure, plot(RMSE0), hold on,
plot(RMSE1),
plot(RMSE2),
% plot(RMSE3),
% plot(RMSE4),
plot(RMSE5),
legend('EnKF','DExKF:IFAC','DExKF:CIN','DExKF:ICI')
% legend('EnKF','DExKF:IFAC','DExKF:CIN','DExKF:CI','DExKF:EI','DExKF:ICI')
x = [1,2,3,4];
avg_time= [avg_time0,avg_time1,avg_time2,avg_time5];
avg_RMSE_cline = [avg_RMSE_cline0,avg_RMSE_cline1,avg_RMSE_cline2,avg_RMSE_cline5];
figure, 
yyaxis left, plot(x,avg_time,'-.*')
ylabel('Average time (sec)')
yyaxis right,plot(x,avg_RMSE_cline,':o')
% axis([D-(D/10) 4*D+(D/10), 0.3 0.5])
ylabel('Average RMSE Cline (sec)');
xlabel('Type of Filter')
xticks(x)
xticklabels({'EnKF','DExKF: IFAC','DExKF: CI(Naive)','DExKF: ICI'})
title('Effect of different fusion algorithms')

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
avg_RMSE_cline0 = RMSE/(Wp.sim.NN - 1)
avg_maxError0 = maxError/(Wp.sim.NN - 1)
avg_time0 = time/(Wp.sim.NN - 1)
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
avg_RMSE_cline1 = RMSE/(Wp.sim.NN - 1)
avg_maxError1 = maxError/(Wp.sim.NN - 1)
avg_time1 = time/(Wp.sim.NN - 1)
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
avg_RMSE_cline2 = RMSE/(Wp.sim.NN - 1)
avg_maxError2 = maxError/(Wp.sim.NN - 1)
avg_time2 = time/(Wp.sim.NN - 1)
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
avg_RMSE_cline3 = RMSE/(Wp.sim.NN - 1)
avg_maxError3 = maxError/(Wp.sim.NN - 1)
avg_time3 = time/(Wp.sim.NN - 1)
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
avg_RMSE_cline4 = RMSE/(Wp.sim.NN - 1)
avg_maxError4 = maxError/(Wp.sim.NN - 1)
avg_time4 = time/(Wp.sim.NN - 1)
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
avg_RMSE_cline5 = RMSE/(Wp.sim.NN - 1)
avg_maxError5 = maxError/(Wp.sim.NN - 1)
avg_time5 = time/(Wp.sim.NN - 1)

figure, plot(RMSE0), hold on,
plot(RMSE1),
plot(RMSE2),
plot(RMSE3),
plot(RMSE4),
plot(RMSE5),
legend('EnKF','DExKF:10','DExKF:20','DExKF:50','DExKF:100','DExKF:Inf')
x = [1,2,3,4,5,6];
avg_time= [avg_time0,avg_time1,avg_time2,avg_time3,avg_time4,avg_time5];
avg_RMSE_cline = [avg_RMSE_cline0,avg_RMSE_cline1,avg_RMSE_cline2,avg_RMSE_cline3,avg_RMSE_cline4,avg_RMSE_cline5];
figure, 
yyaxis left, plot(x,avg_time,'-.*')
ylabel('Average time (sec)')
yyaxis right,plot(x,avg_RMSE_cline,':o')
% axis([D-(D/10) 4*D+(D/10), 0.3 0.5])
ylabel('Average RMSE Cline (sec)');
xlabel('Type of Filter')
xticks(x)
xticklabels({'EnKF','DExKF:10','DExKF:20','DExKF:50','DExKF:100','DExKF:Inf'})
title('Effect of linearizing the system at different time-steps')

%% Disturbances added to the initial flow field (EnKF, DExKF_1D:IFAC,ICI)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_1DZ_noise_init_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim3,timeWFSim3,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    RMSE3(i) = sol_array(i).score.RMSE_cline;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_2DZ_noise_init_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim4,timeWFSim4,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    RMSE4(i) = sol_array(i).score.RMSE_cline;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_dexkf_IFAC_3DZ_noise_init_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim5,timeWFSim5,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    RMSE5(i) = sol_array(i).score.RMSE_cline;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear scriptOptions sol_array strucObs sys Wp
load('C:\Users\Nivas Temp\Documents\Nivas\MSc Thesis\WFObs_Queue\2Turbine\axi_2turb_alm_turb_sim_noise_init_1\workspace.mat')
scriptOptions.powerForecast = 0;
[pwWFSim6,timeWFSim6,~,~] = plotpower( Wp,sol_array,sys,scriptOptions,strucObs, 0 );
for i = 1:Wp.sim.NN
    RMSE6(i) = sol_array(i).score.RMSE_cline;
end
figure, plot(RMSE1), hold on, plot(RMSE2), plot(RMSE3), plot(RMSE4), plot(RMSE5), plot(RMSE6)
legend('EnKF','DExKF_IFAC:1D','DExKF: ICI','DExKF_IFAC:2D','DExKF_IFAC:3D','sim')

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
title('Estimation of Freestream Velcoity'), xlabel('Time (sec)'), ylabel('U_Inf')
figure, plot(RMSE1), hold on, plot(RMSE2)
legend('EnKF','DExKF')
title('RMSE'), xlabel('Time (sec)'), ylabel('RMSE')
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