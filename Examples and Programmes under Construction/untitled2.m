%% Qu = 0.1, Qv = 0.1, R = 0.1
% K: Static, Dynamic 
% U_Inf = 11m/s and it is not estimated 
% (ExKF; EnKF; DExKF: 2D)

clear all
close all
clc

%ExKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf11_noest_Qu0p1Qv0p1R0p1_NLInf/workspace.mat')
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_ExKF_2D_uinf11_noest_Qu0p1Qv0p1R0p1_Pk1k1_HS/workspace.mat')
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf11_noest_Qu0p1Qv0p1R0p1_NL1/workspace.mat')
sol_arrayExKF = sol_array;
xExKF = cell(2000,1);
for i = 1:Wp.sim.NN
    xExKF{i} = sol_arrayExKF(i).x(strucObs.obs_array);
end

% %EnKF: 1D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear scriptOptions sol_array strucObs sys Wp
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_EnKF_1D_uinf11_noest_Qu0p1Qv0p1R0p1_en450_HS/workspace.mat')
% % load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_EnKF_off_uinf11_noest_Qu0p1Qv0p1R0p1_en450/workspace.mat')
% sol_arrayEnKF = sol_array;
% xEnKF = cell(2000,1);
% for i = 1:Wp.sim.NN
%     xEnKF{i} = sol_arrayEnKF(i).x(strucObs.obs_array);
% end
% 
% % %DExKF: 2D: Static%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % clear scriptOptions sol_array strucObs sys Wp
% % load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_2D_uinf11_noest_CIN_Qu0p1Qv0p1R0p1_P20_Static/workspace.mat')
% % sol_arrayDExKFs = sol_array;
% % sol_array2DExKFs = sol_array2;
% % 
% % xDExKFs = cell(2000,1);
% % for i = 1:Wp.sim.NN
% %     xDExKFs{i} = sol_arrayDExKFs(i).x(strucObs.obs_array);
% % end

%DExKF: 2D: Dynamic%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_2D_uinf11_noest_CIN_Qu0p1Qv0p1R0p1_P20_Dyn/workspace.mat')
% load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_2D_uinf11_noest_ICI0p5_Qu0p1Qv0p1R0p1_P5/workspace.mat')
sol_arrayDExKFd = sol_array;
sol_array2DExKFd = sol_array2;

xDExKFd = cell(2000,1);
for i = 1:Wp.sim.NN
    xDExKFd{i} = sol_arrayDExKFd(i).x(strucObs.obs_array);
end

%DExKF2: Fusion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/turb2axialm_ExKF_uinf11_noest_ICI0p5_turbC_Qu0p1Qv0p1R0p1/workspace.mat')
sol_arrayDExKF2 = d_sol_array; 
strucObsDExKF2  = d_strucObs; 

%DExKF3: Fusion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clear scriptOptions sol_array strucObs sys Wp
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_fullmesh_fullest_Queue/turb2axialm_DExKF_2D_uinf11_noest_ICI0p5_turbC_Qu0p1Qv0p1R0p1/workspace.mat')
sol_arrayDExKF3 = d_sol_array; 
strucObsDExKF3  = d_strucObs; 

%DExKF4: No Fusion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/dWFObs_Queue/turb2axialm_4D_ExKF_nofus_turb2_uinf11_noest_Qu0p1Qv0p1R0p1/workspace.mat')
sol_arrayDExKF4 = d_sol_array; 
strucObsDExKF4  = d_strucObs; 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = cell(2000,1);
for i = 1:Wp.sim.NN
    y{i} = sol_array2(i).measuredData.sol(strucObs.obs_array);
end
% [7,11,13,16]
for output = 11
% output  = 15
if output <=5
    k = 1;
    output2 = output
elseif (output>=6)&&(output<=10)
    k = 2;
    output2 = output - 5;
elseif (output>=11)&&(output<=14)
    k = 1;
    output2 = output - 5;
else
    k = 2;
    output2 = output - 9;
end

xDExKF2 = cell(2000,1);
for i = 1:Wp.sim.NN
    xDExKF2{i} = sol_arrayDExKF2{k}(i).x(strucObsDExKF2{k}.obs_array);
end
xDExKF3 = cell(2000,1);
for i = 1:Wp.sim.NN
    xDExKF3{i} = sol_arrayDExKF3{k}(i).x(strucObsDExKF3{k}.obs_array);
end
xDExKF4 = cell(2000,1);
for i = 1:Wp.sim.NN
    xDExKF4{i} = sol_arrayDExKF4{k}(i).x(strucObsDExKF4{k}.obs_array);
end
    
for i = 1:Wp.sim.NN
    yy(i) = y{i}(output);
end
yy=yy';
for i = 1:Wp.sim.NN
    xxExKF(i) = xExKF{i}(output);
end
% for i = 1:Wp.sim.NN
%     xxEnKF(i) = xEnKF{i}(output);
% end
% % for i = 1:Wp.sim.NN
% %     xxDExKFs(i) = xDExKFs{i}(output);
% % end
for i = 1:Wp.sim.NN
    xxDExKFd(i) = xDExKFd{i}(output);
end
for i = 1:Wp.sim.NN
    xxDExKF2(i) = xDExKF2{i}(output2);
end
for i = 1:Wp.sim.NN
    xxDExKF3(i) = xDExKF3{i}(output2);
end
for i = 1:Wp.sim.NN
    xxDExKF4(i) = xDExKF4{i}(output2);
end

% ExKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xxExKF=xxExKF';
XCExKF=xcorr(yy-xxExKF);
% % EnKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% xxEnKF=xxEnKF';
% XCEnKF=xcorr(yy-xxEnKF);
% % % DExKFs%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % xxDExKFs=xxDExKFs';
% % XCDExKFs=xcorr(yy-xxDExKFs);
% DExKFd%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xxDExKFd=xxDExKFd';
XCDExKFd=xcorr(yy-xxDExKFd);
% DExKF2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xxDExKF2=xxDExKF2';
XCDExKF2=xcorr(yy-xxDExKF2);
% DExKF3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xxDExKF3=xxDExKF3';
XCDExKF3=xcorr(yy-xxDExKF3);
% DExKF4%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xxDExKF4=xxDExKF4';
XCDExKF4=xcorr(yy-xxDExKF4);

figure, 
% subplot(2,1,1)
plot(-floor(length(XCExKF)/2):floor(length(XCExKF)/2),XCExKF,'LineWidth',2), 
hold on,
% plot(-floor(length(XCEnKF)/2):floor(length(XCEnKF)/2),XCEnKF,'LineWidth',2)
% % plot(-floor(length(XCDExKFs)/2):floor(length(XCDExKFs)/2),XCDExKFs,'LineWidth',2)
plot(-floor(length(XCDExKFd)/2):floor(length(XCDExKFd)/2),XCDExKFd,'LineWidth',2)
plot(-floor(length(XCDExKF2)/2):floor(length(XCDExKF2)/2),XCDExKF2,'LineWidth',2)
plot(-floor(length(XCDExKF3)/2):floor(length(XCDExKF3)/2),XCDExKF3,'LineWidth',2)
plot(-floor(length(XCDExKF4)/2):floor(length(XCDExKF4)/2),XCDExKF4,'LineWidth',2)
% legend('ExKF','EnKF: 1D','DExKF1_static: 2D','DExKF1_dynamic: 2D','DExKF2','DExKF3: 2D','DExKF4')
legend('Centralized Arch.','Distributed Arch. 1',...
    'Distributed Arch. 2','Distributed Arch. 3','Distributed Arch. 4')
% legend('ExKF','EnKF: 1D','DExKF1: 2D','DExKF2','DExKF3: 2D','DExKF4')
xlabel('Lag (\tau)')
title('Auto correlation of the innovation signal');

% [Xps_ExKF,omExKF]=pwelch(yy-xxExKF,128,[],[],1/1);
% % [Xps_EnKF,omEnKF]=pwelch(yy-xxEnKF,128,[],[],1/1);
% % % [Xps_DExKFs,omDExKFs]=pwelch(yy-xxDExKFs,128,[],[],1/1);
% [Xps_DExKFd,omDExKFd]=pwelch(yy-xxDExKFd,128,[],[],1/1);
% [Xps_DExKF2,omDExKF2]=pwelch(yy-xxDExKF2,128,[],[],1/1);
% [Xps_DExKF3,omDExKF3]=pwelch(yy-xxDExKF3,128,[],[],1/1);
% [Xps_DExKF4,omDExKF4]=pwelch(yy-xxDExKF4,128,[],[],1/1);
% 
% figure,
% % subplot(2,1,2), 
% semilogy(omExKF,Xps_ExKF,'LineWidth',2);
% hold on,
% % semilogy(omEnKF,Xps_EnKF,'LineWidth',2);
% % % semilogy(omDExKFs,Xps_DExKFs,'LineWidth',2);
% semilogy(omDExKFd,Xps_DExKFd,'LineWidth',2);
% semilogy(omDExKF2,Xps_DExKF2,'LineWidth',2);
% semilogy(omDExKF3,Xps_DExKF3,'LineWidth',2);
% semilogy(omDExKF4,Xps_DExKF4,'LineWidth',2);
% % legend('ExKF','EnKF: 1D','DExKF1_static: 2D','DExKF1_dynamic: 2D','DExKF2','DExKF3: 2D','DExKF4')
% legend('Centralized Arch.','Distributed Arch. 1',...
%     'Distributed Arch. 2','Distributed Arch. 3','Distributed Arch. 4')
% % legend('ExKF','EnKF: 1D','DExKF1: 2D','DExKF2','DExKF3: 2D','DExKF4')
% xlabel('\omega (Hz)')
% title('Power Spectrum')
% set(gca,'FontSize')
end