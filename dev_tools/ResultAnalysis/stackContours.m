clear all; clc;

% Time settings
% k_now_vec  = [10 150  300 750]; % Starting point of forecast
k_now_vec  = [1997]; % Starting point of forecast

% Data sources
data = {};

% COMPARE ExKF: NL1, NL10, NL20, NL50, NL100, NLInf
% data{end+1} = struct(...
%     'name','Frequency: 1 Hz',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf11_noest_Qu0p1Qv0p1R0p1_NL1/workspace.mat');
% data{end+1} = struct(...
%     'name','Frequency: 1/10 Hz',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf11_noest_Qu0p1Qv0p1R0p1_NL10/workspace.mat');
% data{end+1} = struct(...
%     'name','Frequency: 1/20 Hz',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf11_noest_Qu0p1Qv0p1R0p1_NL20/workspace.mat');
% data{end+1} = struct(...
%     'name','Frequency: 1/50 Hz',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf11_noest_Qu0p1Qv0p1R0p1_NL50/workspace.mat');
% data{end+1} = struct(...
%     'name','Frequency: 1/100 Hz',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf11_noest_Qu0p1Qv0p1R0p1_NL100/workspace.mat');
% data{end+1} = struct(...
%     'name','Only at the first iteration',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf11_noest_Qu0p1Qv0p1R0p1_NLInf/workspace.mat');
% outputFigName = ['k' strrep(num2str(k_now_vec),' ','_') '_flowContours_SOWFAComparison'];

% data{end+1} = struct(...
%     'name','Localization size: 1D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_ExKF_1D_uinf11_noest_Qu0p1Qv0p1R0p1_Pk1k1_HS/workspace.mat');
% data{end+1} = struct(...
%     'name','Localization size: 2D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_ExKF_2D_uinf11_noest_Qu0p1Qv0p1R0p1_Pk1k1_HS/workspace.mat');
% data{end+1} = struct(...
%     'name','Localization size: 3D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_ExKF_3D_uinf11_noest_Qu0p1Qv0p1R0p1_Pk1k1_HS/workspace.mat');
% data{end+1} = struct(...
%     'name','Localization size: 4D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_ExKF_4D_uinf11_noest_Qu0p1Qv0p1R0p1_Pk1k1_HS/workspace.mat');
% data{end+1} = struct(...
%     'name','Localization size: Full',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/axi_2turb_alm_turb_ExKF_uinf11_noest_Qu0p1Qv0p1R0p1_NLInf/workspace.mat');
% outputFigName = ['k' strrep(num2str(k_now_vec),' ','_') '_flowContours_SOWFAComparison'];

% EnKF
% data{end+1} = struct(...
%     'name','Localization size: 1D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_EnKF_1D_uinf11_noest_Qu0p1Qv0p1R0p1_en150_HS/workspace.mat');
% data{end+1} = struct(...
%     'name','Localization size: 2D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_EnKF_2D_uinf11_noest_Qu0p1Qv0p1R0p1_en450_HS/workspace.mat');
% data{end+1} = struct(...
%     'name','Localization size: 3D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_EnKF_3D_uinf11_noest_Qu0p1Qv0p1R0p1_en450_HS/workspace.mat');
% data{end+1} = struct(...
%     'name','Localization size: 4D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_EnKF_4D_uinf11_noest_Qu0p1Qv0p1R0p1_en450_HS/workspace.mat');
% data{end+1} = struct(...
%     'name','Localization size: Full',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_EnKF_off_uinf11_noest_Qu0p1Qv0p1R0p1_en450/workspace.mat');
% outputFigName = ['k' strrep(num2str(k_now_vec),' ','_') '_flowContours_SOWFAComparison'];

%DExKF
data{end+1} = struct(...
    'name','Domain size: 1D',...
    'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_1D_uinf11_noest_ICI0p5_Qu0p1Qv0p1R0p1_P20/workspace.mat');
data{end+1} = struct(...
    'name','Domain size: 2D',...
    'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_2D_uinf11_noest_ICI0p5_Qu0p1Qv0p1R0p1_P20/workspace.mat');
data{end+1} = struct(...
    'name','Domain size: 3D',...
    'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_3D_uinf11_noest_ICI0p5_Qu0p1Qv0p1R0p1_P20/workspace.mat');
% data{end+1} = struct(...
%     'name','Domain size: 4D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_4D_uinf11_noest_ICI0p5_Qu0p1Qv0p1R0p1_P20/workspace.mat');
outputFigName = ['k' strrep(num2str(k_now_vec),' ','_') '_flowContours_SOWFAComparison'];

%DExKF: NO fusion
% data{end+1} = struct(...
%     'name','Domain size: 1D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_1D_uinf11_noest_no_Qu0p1Qv0p1R0p1_P20/workspace.mat');
% data{end+1} = struct(...
%     'name','Domain size: 2D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_2D_uinf11_noest_no_Qu0p1Qv0p1R0p1_P20/workspace.mat');
% data{end+1} = struct(...
%     'name','Domain size: 3D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_3D_uinf11_noest_no_Qu0p1Qv0p1R0p1_P20/workspace.mat');
% data{end+1} = struct(...
%     'name','Domain size: 4D',...
%     'path','/Users/Nivas_Kumar/Documents/NivasStudyMaterials/TUDelft/EnKF+WFSim/WFObs_Queue/2Turbine/turb2axialm_DExKF_4D_uinf11_noest_no_Qu0p1Qv0p1R0p1_P20/workspace.mat');
% outputFigName = ['k' strrep(num2str(k_now_vec),' ','_') '_flowContours_SOWFAComparison'];
%% Core operations
addpath('../../WFSim/libraries/export_fig'); % Add export_fig library
addpath('../../dev_tools/ResultAnalysis/libraries'); % Add libraries for VAF calculations

for di = 1:length(data)
    % Load workspace file
    disp(['' num2str(di) '. Loading workspace.mat for ''' data{di}.name '''.']);
    WS{di} = load(data{di}.path);
    
    for kn = 1:length(k_now_vec)
        k_now = k_now_vec(kn);
        sol   = WS{di}.sol_array(k_now);
        
        out(kn,di).u  = sol.u;
        out(kn,di).uq = sol.measuredData.uq;
        out(kn,di).e  = abs(sol.u-sol.measuredData.uq);
    end
end

%% Plot figures
% Meshing settings
x = WS{1}.Wp.mesh.ldyy;
y = WS{1}.Wp.mesh.ldxx2;

% Plot velocities in a contourf figure
close all; h = figure; 
% if length(data) == 3
%     h.Position = [448.2000 52.2000 778.4000 726.4000];
% end
% if length(data) == 4
%     h.Position = [1.5666e+03 -19.8000 802.4000 967.2000];
% end
% set(h,'defaultTextInterpreter','latex')
% set(0,'defaultTextInterpreter','latex')
climits_u = [0 12];
climits_e = [0 3];
nF_vert = length(data); % number of figures vertically
nF_horz = 3;  % number of figures horizontally

% applied correction for yaw angle: wake was forming at wrong side
for j = 1:length(data)
	rotorRotation = -.5*WS{1}.Wp.turbine.Drotor*exp(1i*-WS{1}.Wp.turbine.input(k_now_vec).phi'*pi/180);
        subaxis(nF_vert,nF_horz,3*(j-1)+1,'SpacingVert',0.01,'SpacingHoriz',0.00);
        contourf(x,y,out(1,j).uq,[0:0.1:12],'Linecolor','none');
        caxis(climits_u);
        if j ~= length(data)
            set(gca,'XTickLabel',[]);
        end
        set(gca,'YTick',0:500:max(y(:)));
        ylabel({data{j}.name;'x-direction (m)'})
        if j == length(data)
            xlabel('y-dir. (m)')
        end
        axis equal tight;
        if j == 1
            title('SOWFA');
        end
        turb(WS,rotorRotation);
        subaxis(nF_vert,nF_horz,3*(j-1)+2,'SpacingVert',0.01,'SpacingHoriz',0.00);
        contourf(x,y,out(1,j).u,[0:0.1:12],'Linecolor','none');
        caxis(climits_u);
        if j ~= length(data)
            set(gca,'XTickLabel',[]);
        end
        set(gca,'YTickLabel',[]);
        axis equal tight;
        if j == 1
%             title('EnKF');
            title('DExKF: Architecture 1');
        end
        if j == length(data)
            xlabel('y-dir. (m)')
        end
        turb(WS,rotorRotation);
        if j == length(data)
            xlabel('y-dir. (m)')
            clb_u = colorbar('south');
            clb_u.Position = [0.1625 0.0475 0.4000 0.0100];
            clb_u.Label.String = 'Flow speed (m/s)';
            caxis(climits_u);
        end        
        subaxis(nF_vert,nF_horz,3*(j-1)+3,'SpacingVert',0.01,'SpacingHoriz',0.00);
        contourf(x,y,out(1,j).e,[0:0.1:3],'Linecolor','none');
        caxis(climits_e);
        if j ~= length(data)
            set(gca,'XTickLabel',[]);
        end
        set(gca,'YTickLabel',[]);
        axis equal tight;
        if j == 1
            title('Error');
        end
        if j == length(data)
%             xlabel('y-dir. (m)')
%             clb_u = colorbar('south');
%             clb_u.Position = [0.1625 0.0475 0.4000 0.0100];
%             clb_u.Label.String = 'Flow speed (m/s)';
%             caxis(climits_u);
            clb_e = colorbar('south');
            clb_e.Position = [0.6700 0.0475 0.2300 0.0100];
            clb_e.Label.String = 'Error (m/s)';
            clb_e.Limits = climits_e;
        end        

        % Turbines
        turb(WS,rotorRotation);
end
colormap(jet)


% export_fig(outputFigName,'-pdf','-transparent')


function turb(WS,rotorRotation)
        hold all;
        Wp = WS{1}.Wp;
        for kk=1:Wp.turbine.N
            Qy     = (Wp.turbine.Cry(kk)-abs(real(rotorRotation(kk)))):1:(Wp.turbine.Cry(kk)+abs(real(rotorRotation(kk))));
            Qx     = linspace(Wp.turbine.Crx(kk)-imag(rotorRotation(kk)),Wp.turbine.Crx(kk)+imag(rotorRotation(kk)),length(Qy));
            rectangle('Position',[Wp.turbine.Cry(kk)-0.10*Wp.turbine.Drotor Wp.turbine.Crx(kk) ...
                0.20*Wp.turbine.Drotor 0.30*Wp.turbine.Drotor],'Curvature',0.2,...
                'FaceColor','w')
            plot(mean(Qy)+(Qy-mean(Qy))*1.2,mean(Qx)+(Qx-mean(Qx))*1.2,'k','linewidth',3)
            plot(Qy,Qx,'w','linewidth',2)
        end
            % Sensors
%             if WS{j}.strucObs.measFlow
%                 plot([WS{j}.strucObs.obs_array_locu.y],[WS{j}.strucObs.obs_array_locu.x],'w.','lineWidth',3.0,'displayName','Sensors');
%                 plot([WS{j}.strucObs.obs_array_locu.y],[WS{j}.strucObs.obs_array_locu.x],'r.','displayName','Sensors');
%             end
        set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix

end