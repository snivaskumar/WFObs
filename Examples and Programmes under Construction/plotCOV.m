function plotCOV( Wp,sol_array,sys,scriptOptions,strucObs);

sol                     = sol_array(end);
Nx                      = Wp.mesh.Nx;
Ny                      = Wp.mesh.Ny;
if ~length(strfind(scriptOptions.savePath,'EnKF'))
    P       = diag(strucObs.Pk);
else
    Aenf    = strucObs.Aen;
    Aenft   = Aenf-repmat(mean(Aenf,2),1,strucObs.nrens);
%     P       = ( 1/(strucObs.nrens-1) )*Aenft*Aenft';
    P       = Aenft*Aenft';
    P       = diag(P);
end
tmp                     = max(max(P));
Puv                     = tmp*ones(Nx,Ny);
Puv(3:end-1,2:end-1)    = reshape(P(1:(Nx-3)*(Ny-2)),Ny-2,Nx-3)';
    
hFigs = {};
scrsz = get(0,'ScreenSize'); 
hFigs{1}=figure('color',[1 1 1],'Position',[50 50 floor(scrsz(3)/1.1) floor(scrsz(4)/1.1)], 'MenuBar','none','ToolBar','none','visible', 'on');
set(hFigs{1},'defaultTextInterpreter','latex')
data{1} = struct('x',Wp.mesh.ldyy, 'y',Wp.mesh.ldxx2,'z',Puv,'title','Co_variance');
        
% applied correction for yaw angle: wake was forming at wrong side
rotorRotation = -.5*Wp.turbine.Drotor*exp(1i*-Wp.turbine.input(sol.k).phi'*pi/180); 
    
% Plot velocities in a contourf figure
set(0,'CurrentFigure',hFigs{1}); clf
subplotDim = numSubplots(length(data));
for j = 1:length(data)
	subplot(subplotDim(1),subplotDim(2),j);
    cmax = (max(max(Puv)));
    cmin = (min(min(Puv)));
    contourf(data{j}.x,data{j}.y,data{j}.z,cmin:0.1:cmax,'Linecolor','none');
    title([data{j}.title ' (t = ' num2str(sol.time) ')'])
    hold all; colorbar;
    caxis([cmin cmax]);
    axis equal; axis tight;
    xlabel('y-direction')
    ylabel('x-direction')         
    hold all
    % Turbines
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
    if strucObs.measFlow
    	plot([strucObs.obs_array_locu.y],[strucObs.obs_array_locu.x],'wo','lineWidth',3.0,'displayName','Sensors');
        plot([strucObs.obs_array_locu.y],[strucObs.obs_array_locu.x],'ro','displayName','Sensors');
    end
	set(gca,'YDir','Reverse'); % Flip axis so plot matches matrix
end;
colormap(jet)
end