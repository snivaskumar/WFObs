function [Wp,sol_out,strucObs] = WFObs_o_exkf(strucObs,Wp,sys_in,sol_in,options)
% WFOBS_O_EXKF  Extended KF algorithm for recursive state estimation
%
%   SUMMARY
%    This code performs state estimation using the Extended Kalman filter
%    (ExKF) algorithm. It uses high-fidelity measurements
%    (sol.measuredData) to improve the flow estimation compared to
%    open-loop simulations with WFSim.
%
%   RELEVANT INPUT/OUTPUT VARIABLES
%      see 'WFObs_o.m' for the complete list.
%   

% Import measurement data variable
measuredData = sol_in.measuredData;

if sol_in.k == 1
    % Setup covariance and system output matrices
    if options.exportPressures
        strucObs.Pk    = sparse(eye(strucObs.size_state))*strucObs.P_0;
        strucObs.Htt   = sparse(eye(strucObs.size_state));
        strucObs.Htt   = strucObs.Htt(strucObs.obs_array,:);
%         strucObs.Q_k   = strucObs.Q_k*eye(strucObs.size_state);
        Qu             = strucObs.Q_e.u*ones(1,Wp.Nu);
        Qv             = strucObs.Q_e.v*ones(1,Wp.Nv);
        Qp             = strucObs.Q_e.p*ones(1,Wp.Np);
        strucObs.Q_k   = [Qu, Qv, Qp].*eye(strucObs.size_state);
    else
        strucObs.Pk    = sparse(eye(strucObs.size_output))*strucObs.P_0;
        strucObs.Htt   = sparse(eye(strucObs.size_output));
        strucObs.Htt   = strucObs.Htt(strucObs.obs_array,:);
%         strucObs.Q_k   = strucObs.Q_k*eye(strucObs.size_output);
        Qu             = strucObs.Q_e.u*ones(1,Wp.Nu);
        Qv             = strucObs.Q_e.v*ones(1,Wp.Nv);
%         Qp             = strucObs.Q_e.p*ones(1,Wp.Np);
        strucObs.Q_k   = [Qu, Qv].*eye(strucObs.size_output);
    end;
end;

% ExKF forecast update
soltemp   = sol_in;
soltemp.k = soltemp.k - 1;
[solf,sysf]             = WFSim_timestepping( soltemp, sys_in, Wp, options );       % Forward propagation
% Fk(sysf.pRCM,sysf.pRCM) = sysf.A(sysf.pRCM,sysf.pRCM)\sysf.Al(sysf.pRCM,sysf.pRCM); % Linearized A-matrix at time k
% Bk(sysf.pRCM,:)         = sysf.A(sysf.pRCM,sysf.pRCM)\sysf.Bl(sysf.pRCM,:);         % Linearized B-matrix at time k

NL = strucObs.linearize_freq;
if (sol_in.k == 1) || (rem(sol_in.k,NL) == 0)
    clear Fk Bk
    Fk(sysf.pRCM,sysf.pRCM) = sysf.A(sysf.pRCM,sysf.pRCM)\sysf.Al(sysf.pRCM,sysf.pRCM); % Linearized A-matrix at time k
    Bk(sysf.pRCM,:)         = sysf.A(sysf.pRCM,sysf.pRCM)\sysf.Bl(sysf.pRCM,:);         % Linearized B-matrix at time k
    strucObs.Fk = Fk;
    strucObs.Bk = Bk;
end
Fk = strucObs.Fk;
Bk = strucObs.Bk;
n  = length(Fk);

% Neglect pressure terms
if ~options.exportPressures 
    Fk     = Fk(1:strucObs.size_output,1:strucObs.size_output);
    Bk     = Bk(1:strucObs.size_output,:);
    solf.x = solf.x(1:strucObs.size_output);
end;
Pk = strucObs.Pk;
% Pf = Fk*Pk*Fk' + strucObs.Q_k;  % Covariance matrix P for x(k) knowing y(k-1)

if strucObs.localize
%     if sol_in.k == 1
%         if strucObs.stateEst || strucObs.measFlow
%             stateLocArray = zeros(strucObs.size_output,2);
%             for iii = 1:strucObs.size_output
%                 [~,loci,~]           = WFObs_s_sensors_nr2grid(iii,Wp.mesh);
%                 stateLocArray(iii,:) = [loci.x, loci.y];
%             end
%         end
% 
%         % Generate the locations of all turbines
%         if strucObs.tune.est || strucObs.measPw
%             turbLocArray = zeros(Wp.turbine.N,2);
%             for iii = 1:Wp.turbine.N
%                 turbLocArray(iii,:) = [Wp.turbine.Crx(iii),Wp.turbine.Cry(iii)];
%             end
%         end
% 
%         % Generate the locations of all outputs
%         outputLocArray = [];
%         if strucObs.measFlow
%             outputLocArray = [outputLocArray; stateLocArray(strucObs.obs_array,:)];
%         end
%         if strucObs.measPw
%             outputLocArray = [outputLocArray; turbLocArray];
%         end
%         if strucObs.stateEst
%             rho_locl.cross = sparse(strucObs.size_output,strucObs.size_output);
%             for iii = 1:strucObs.size_output % Loop over all default states
%                 loc1 = stateLocArray(iii,:);
%                 for jjj = 1:length(strucObs.obs_array) % Loop over all measurements
%                     loc2 = outputLocArray(jjj,:);
%                     dx = sqrt(sum((loc1-loc2).^2)); % displacement between state and output
%                     if dx <= strucObs.l_locl
%                         strucObs.rho_locl(iii,jjj) = 1;
%                     else
%                         strucObs.rho_locl(iii,jjj) = 0;
%                     end;
%                 end
%             end
%             clear iii jjj dx loc1 loc2
%         else
%             rho_locl.cross = [];
%         end
%     end
%     factor = strucObs.rho_locl;

    strucObs.nrobs    = length(strucObs.obs_array); % number of state measurements
    strucObs.M        = strucObs.nrobs+strucObs.measPw*Wp.turbine.N; % total length of measurements
    if sol_in.k == 1
%         options.savePath
        strucObs.l_locl
        [ strucObs ] = WFObs_o_exkf_localization( Wp,strucObs );
%         options.savePath
        save([options.savePath '/workspace0.mat'],'strucObs');
    end

% Pk = strucObs.autoState_corrfactor.* Pk;
end
Pf = Fk*Pk*Fk' + strucObs.Q_k;  % Covariance matrix P for x(k) knowing y(k-1)
    
% ExKF analysis update
sol_out     = sol_in; % Copy previous solution before updating x

if strucObs.localize && (strucObs.localizeType == 1)
    Kgain       = strucObs.cross_corrfactor.* Pf(:,strucObs.obs_array)...
                  *pinv(strucObs.auto_corrfactor .* Pf(strucObs.obs_array,...
                  strucObs.obs_array)+strucObs.R_k); % Kalman gain           
else
    Kgain       = Pf(:,strucObs.obs_array)*pinv(Pf(strucObs.obs_array,...
                  strucObs.obs_array)+strucObs.R_k); % Kalman gain           
end

sol_out.x   = solf.x + Kgain*(measuredData.sol(strucObs.obs_array)...
                 -solf.x(strucObs.obs_array)); % Optimally predicted state vector
             
if strucObs.localize && (strucObs.localizeType == 2)
    Pf          = strucObs.autoState_corrfactor.* Pf; 
end

strucObs.Pk = (eye(size(Pf))-Kgain*strucObs.Htt)*Pf;  % State covariance matrix

% Pk1k1   = strucObs.Pk;
% R       = strucObs.R_k;
% Q       = strucObs.Q_k; 
% y       = measuredData.sol(strucObs.obs_array); 
% C       = strucObs.Htt; 
% n       = length(Fk);
% 
% xkk1    = solf.x;
% Pkk1    = Fk*Pk1k1*Fk' + Q;
% dy      = y - C*xkk1;
% Pyy     = R + C*Pkk1*C';
% Pxy     = Pkk1*C';
% K       = Pxy*inv(Pyy);
% xkk     = xkk1 + K*dy;
% Pkk     = ( eye(n,n) - K*C )*Pkk1;
% 
% sol_out.x   = xkk;
% strucObs.Pk = Pkk;

% Export new solution from estimation
[sol_out,~]  = MapSolution(Wp,sol_out,Inf,options); % Map solution to flowfields
[~,sol_out]  = Actuator(Wp,sol_out,options);        % Recalculate power after analysis update
end