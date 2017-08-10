clear all; clc; close all;

%% Set-up
% Source files
scriptOptions.plotMapping = false; % Plot mapping at time k == 1  (validation)
scriptOptions.sourcePath  = ['D:/bmdoekemeijer/My Documents/MATLAB/WFObs/WFSim/data_PALM/2turb_adm_long'];

% Turbine properties directly from PALM
rawTurbData.Crx       = [5700, 6456];   % Raw (inertial frame) in (m)
rawTurbData.Cry       = [1175, 1175];   % Raw (inertial frame) in (m)
rawTurbData.hubHeight = 90.0;           % Hub height in (m)

% Desired output settings
meshSetup.name        = 'lin50x25';
meshSetup.dt          = 1.0 ; % Timestep in seconds
meshSetup.distance_S  = 300 ; % distance (m) upwind   first  turbine to export
meshSetup.distance_N  = 740;  % distance (m) downwind  last  turbine to export
meshSetup.distance_W  = 280 ; % distance (m) west most left  turbine to export
meshSetup.distance_E  = 280 ; % distance (m) east most right turbine to export
meshSetup.Nx          = 50;
meshSetup.Ny          = 25;


%% Core code
tic

% Sort and import files from PALM data folder
disp('Sorting and importing files from source folder...')
filesInFolder = dir(scriptOptions.sourcePath); % Gather all files from PALM/SOWFA folder
filesInFolder = {filesInFolder(3:end).name};   % Remove '.' and '..'
for jFile = filesInFolder
    jFile = jFile{1};
    % Check if it is a flow data
    if ~isempty(strfind(jFile,'.nc')) % if contains .nc
        if exist('flowDataRaw')
            error('Confused by the existence of multiple .nc files.');
        else
            flowDataRaw = loadFlowData(jFile);
        end
        % Otherwise, check if it is a turbine file
    elseif ~isempty(strfind(jFile,'turbine_param')) % if contains 'turbine_param'
        stringIdx = strfind(jFile,'turbine_param');
        turbineId = regexp(jFile(stringIdx:end),'\d*','Match');
        turbineId = str2double(turbineId); % Convert to turbine Id

        % Import from external file
        turbDataRaw(turbineId) = loadTurbData(jFile);
    end
end
clear filesInFolder jFile stringIdx turbineId % Clear variables

% Convert array of structs 'turbData' to single struct 'turbData'
fieldNamesListTurb = fieldnames(turbDataRaw);
for j = 1:length(fieldNamesListTurb) 
    jName = fieldNamesListTurb{j};
    turbDataTmp.(jName) = [turbDataRaw.(jName)];
end
turbDataRaw = turbDataTmp;
clear j jName turbDataTmp % Clear variables

% Resample data
disp('Resampling data in time...')
time  = meshSetup.dt:meshSetup.dt:min([flowDataRaw.time(end),turbDataRaw.t(end)]);
k_raw = interp1([0; flowDataRaw.time],[1 1:length(flowDataRaw.time)],time,'nearest');
flowDataResampled      = flowDataRaw;
flowDataResampled.time = time;
for j = {'u','v','w'} % Resample flow fields
    flowDataResampled.(j{1}) = flowDataRaw.(j{1})(k_raw,:,:,:);
end
k_raw = interp1(turbDataRaw.t(:,1),1:length(turbDataRaw.t(:,1)),time,'nearest');
for j = 1:length(fieldNamesListTurb)
    turbDataResampled.(fieldNamesListTurb{j}) = turbDataRaw.(fieldNamesListTurb{j})(k_raw,:);
end
turbDataResampled.t = time;
clear j k_raw


% Normalize variables
flowData     = flowDataResampled;
flowData.x   = flowData.x  - flowData.x(1);
flowData.y   = flowData.y  - flowData.y(1);
flowData.xu  = flowData.xu - flowData.xu(1);
flowData.yv  = flowData.yv - flowData.yv(1);

turbData        = turbDataResampled;
turbData.Crx    = rawTurbData.Crx - flowDataRaw.xu(1);
turbData.Cry    = rawTurbData.Cry - flowDataRaw.yv(1);

% Find z-index closest to turbine hub-height
nz = round(interp1(flowData.zw_3d,1:length(flowData.zw_3d),rawTurbData.hubHeight));

% Determine freestream conditions flow field by checking all corners
disp('Rotating and translating grid according to u_Inf and v_Inf...')
U_Inf = 0;
for i = [1 size(flowData.u,3)]
    for j = [1 size(flowData.u,4)]
        tmp_u = mean(flowData.u(:,nz,i,j));
        tmp_v = mean(flowData.v(:,nz,i,j));
        tmp_U = sqrt(tmp_u^2 + tmp_v^2);
        if U_Inf < tmp_U
            u_Inf = tmp_u; 
            v_Inf = tmp_v;
            U_Inf = tmp_U;
        end
    end
end
clear tmp_u tmp_v tmp_U i j

% Rotate wind field
windDirection = atan(v_Inf/u_Inf); % in radians
if windDirection > deg2rad(2.5)
    % Rotate flow field
    error('Code not yet written -- has not appeared necessary yet.');
    % ...
end

% Preprocessing for translation: check if even possible
if meshSetup.distance_S > min(turbData.Crx)
    disp(' ')
    disp('WARNING: Specified output settings for x exceed available data.')
    disp('Shrinking to a smaller distance_S that matches the available data.')
    meshSetup.distance_S = min(turbData.Crx)
end
yMax  = flowData.y(end)-flowData.y(1);
xMax  = flowData.x(end)-flowData.x(1);
xTurbSeperation = max(turbData.Crx)-min(turbData.Crx);
Wp.Lx = meshSetup.distance_S + meshSetup.distance_N + xTurbSeperation;
Wp.Ly = meshSetup.distance_W + meshSetup.distance_E;
if Wp.Lx > xMax
    disp(' ')
    disp('WARNING: Specified output settings for x exceed available data.')
    disp('Resizing to a symmetric mesh that matches the available data.')
    [meshSetup.distance_S,meshSetup.distance_N] = deal((xMax - xTurbSeperation) / 2)
end
if Wp.Ly > yMax
    disp('WARNING: Specified output settings for y exceed available data.');
    disp('Resizing to a symmetric mesh that matches the available data.');
    [meshSetup.distance_W,meshSetup.distance_E] = deal(yMax / 2)
end
clear xMax yMax xTurbSeperation

% Translation
flowData.x  = flowData.x  - min(turbData.Crx) + meshSetup.distance_S;
flowData.xu = flowData.xu - min(turbData.Crx) + meshSetup.distance_S;
flowData.y  = flowData.y  - min(turbData.Cry) + meshSetup.distance_W;
flowData.yv = flowData.yv - min(turbData.Cry) + meshSetup.distance_W;
        
% Determine target meshing (simplified meshing.m code)
Wp.ldx   = linspace(0,Wp.Lx,meshSetup.Nx);
Wp.ldy   = linspace(0,Wp.Ly,meshSetup.Ny);
Wp.ldxx  = repmat(Wp.ldx',1,meshSetup.Ny);
Wp.ldyy  = repmat(Wp.ldy,meshSetup.Nx,1);
Wp.ldx2  = 0.5*(Wp.ldx(1:end-1)+Wp.ldx(2:end));
Wp.ldx2  = [Wp.ldx2 2*Wp.ldx2(end)-Wp.ldx2(end-1)]; 
Wp.ldy2  = 0.5*(Wp.ldy(1:end-1)+Wp.ldy(2:end));
Wp.ldy2  = [Wp.ldy2 2*Wp.ldy2(end)-Wp.ldy2(end-1)];
Wp.ldxx2 = repmat(Wp.ldx2',1,meshSetup.Ny);
Wp.ldyy2 = repmat(Wp.ldy2,meshSetup.Nx,1);
Wp.Nx    = meshSetup.Nx;
Wp.Ny    = meshSetup.Ny;


%% Perform remesh for every timestep
NN      = length(time);
Nx_raw  = size(flowData.u,3);
Ny_raw  = size(flowData.u,4);
[yu,xu] = meshgrid(flowData.y,flowData.xu);
[yv,xv] = meshgrid(flowData.yv,flowData.x);
[u,v]   = deal(zeros(NN,meshSetup.Nx,meshSetup.Ny));
for k = 1:5%NN
    disp(['Performing remesh for k = ' num2str(k) '...']);
    
    % u-velocity
    uk_raw      = reshape(flowData.u(k,nz,:,:),Nx_raw,Ny_raw)';  % u(k,z,y,x)   
    uk_remeshed = griddata(yu,xu,uk_raw,Wp.ldyy,Wp.ldxx2, 'nearest');
    uk_remeshed(isnan(uk_remeshed)) = 0; % remove NaNs
    
    % v-velocity
    vk_raw      = reshape(flowData.v(k,nz,:,:),Nx_raw,Ny_raw)';  % u(k,z,y,x)   
    vk_remeshed = griddata(yv,xv,vk_raw,Wp.ldyy2,Wp.ldxx, 'nearest');
    vk_remeshed(isnan(uk_remeshed)) = 0; % remove NaNs    
    
    % save to tensors
    u(k,:,:) = uk_remeshed;
    v(k,:,:) = vk_remeshed;
        
    if k == 1 && scriptOptions.plotMapping
        clf;
        subplot(1,2,1);
        contourf(flowDataRaw.y,flowDataRaw.xu,uk_raw,'Linecolor','none');  
        caxis([min(min(uk_raw)) max(max(uk_raw))+.01]);  hold on; 
        axis equal; axis tight; colormap(hot); colorbar;
        for jTurb = 1:length(turbData.Crx)
            plot(rawTurbData.Cry(jTurb)+[-60,60],rawTurbData.Crx(jTurb)*[1,1],'k-');
        end
        ylabel('x-direction (m)'); xlabel('y-direction (m)');
        title('RAW $u$ [m/s]','interpreter','latex');

        subplot(1,2,2);
        contourf(Wp.ldyy,Wp.ldxx2,uk_remeshed,'Linecolor','none');  colormap(hot);
        caxis([min(min(uk_raw)) max(max(uk_raw))+.01]);  hold on; 
        axis equal; axis tight; colormap(hot); colorbar;
        for jTurb = 1:length(turbData.Crx)
            plot(turbData.Cry(jTurb)+[-60,60],turbData.Crx(jTurb)*[1,1],'k-');
        end
        ylabel('x-direction (m)'); xlabel('y-direction (m)');
        title('RESIZED $u$ [m/s]','interpreter','latex');
        drawnow;    
        
        clear jTurb 
    end
end
clear Nx_raw Ny_raw k uk_raw uk_remeshed vk_raw vk_remeshed

% Save output data
disp('Saving output data...');
[fname,pth] = uiputfile('.mat');
save([pth fname],'time','u','v','xu','yu','xv','yv','turbData','Wp')
toc
 
% Functions for loading datafiles
function turbDataOut = loadTurbData(fileName)
M = load(fileName);
% [Time   UR  Uinf  Ct_adm  a Yaw Thrust Power  WFPower]
turbDataOut.t    = M(:,1);
turbDataOut.Ur   = M(:,2);
turbDataOut.Uinf = M(:,3);
turbDataOut.CT   = M(:,4);
turbDataOut.a    = M(:,5);
turbDataOut.phi  = M(:,6);
turbDataOut.P    = M(:,8);
end

function flowDataOut = loadFlowData(fileName)
flowDataOut.time  = double(nc_varget(fileName,'time'));
flowDataOut.u     = double(nc_varget(fileName,'u'));
flowDataOut.v     = double(nc_varget(fileName,'v'));
flowDataOut.w     = double(nc_varget(fileName,'w'));
flowDataOut.x     = double(nc_varget(fileName,'x'));
flowDataOut.y     = double(nc_varget(fileName,'y'));
flowDataOut.xu    = double(nc_varget(fileName,'xu'));
flowDataOut.yv    = double(nc_varget(fileName,'yv'));
flowDataOut.zw_3d = double(nc_varget(fileName,'zw_3d'));
flowDataOut.zu_3d = double(nc_varget(fileName,'zu_3d'));
flowDataOut.dx    = diff(flowDataOut.x);
flowDataOut.dy    = diff(flowDataOut.y);
flowDataOut.dxu   = diff(flowDataOut.xu);
flowDataOut.dyv   = diff(flowDataOut.yv);
end