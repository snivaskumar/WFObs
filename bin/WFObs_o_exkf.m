function [Wp,sol_out,strucObs] = WFObs_o_exkf(strucObs,Wp,sys_in,sol_in,options,sol_array2)
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
    strucObs.nrobs    = length(strucObs.obs_array); % number of state measurements
    strucObs.M        = strucObs.nrobs+strucObs.measPw*Wp.turbine.N; % total length of measurements
    if sol_in.k == 1
        strucObs.l_locl
        [ strucObs ] = WFObs_o_exkf_localization( Wp,strucObs );
        save([options.savePath '/workspace0.mat'],'strucObs');
    end

end
if strucObs.localize && (strucObs.localizeType == 3)
    Pk          = strucObs.autoState_corrfactor.* Pk; 
end

if ( strcmp(upper(strucObs.KF),'DYNAMIC') )...
        ||( sol_in.k<=2 )...
        ||(( sol_array2(sol_in.k-1).sPk - sol_array2(sol_in.k-2).sPk )>=1 )
    Pf = Fk*Pk*Fk' + strucObs.Q_k;  % Covariance matrix P for x(k) knowing y(k-1)
end
    
if strucObs.localize && (strucObs.localizeType == 4)
    Pf          = strucObs.autoState_corrfactor.* Pf; 
end

% ExKF analysis update
sol_out     = sol_in; % Copy previous solution before updating x

if ( strcmp(upper(strucObs.KF),'DYNAMIC') )...
        ||( sol_in.k<=2 )...
        ||(( sol_array2(sol_in.k-1).sPk - sol_array2(sol_in.k-2).sPk )>=1 )
    if strucObs.localize && (strucObs.localizeType == 1)
        Kgain       = strucObs.cross_corrfactor.* Pf(:,strucObs.obs_array)...
                      *pinv(strucObs.auto_corrfactor .* Pf(strucObs.obs_array,...
                      strucObs.obs_array)+strucObs.R_k); % Kalman gain           
    else
        Kgain       = Pf(:,strucObs.obs_array)*pinv(Pf(strucObs.obs_array,...
                      strucObs.obs_array)+strucObs.R_k); % Kalman gain           
    end
else 
    Kgain = sol_array2(sol_in.k-1).K;
end

sol_out.x   = solf.x + Kgain*(measuredData.sol(strucObs.obs_array)...
                 -solf.x(strucObs.obs_array)); % Optimally predicted state vector
             
if strucObs.localize && (strucObs.localizeType == 2)
    Pf          = strucObs.autoState_corrfactor.* Pf; 
end

if ( strcmp(upper(strucObs.KF),'DYNAMIC') )...
        ||( sol_in.k<=2 )...
        ||(( sol_array2(sol_in.k-1).sPk - sol_array2(sol_in.k-2).sPk )>=1 )
    strucObs.Pk = (eye(size(Pf))-Kgain*strucObs.Htt)*Pf;  % State covariance matrix
end

sol_out.K = Kgain;

% tmp1 = (abs(eig(Pf)));
% sol_out.Pkk1 = tmp1;
% % plot(sol_out.k,tmp1,'o'),hold on;
% stmp1 = max(abs(eig(Pf)));
% sol_out.sPkk1 = stmp1;

Pk = strucObs.Pk;
tmp2 = (abs(eig(Pk)));
% plot(sol_out.k,tmp2,'*'),hold on;
sol_out.Pk = tmp2;
stmp2 = max(abs(eig(Pk)));
sol_out.sPk = stmp2;

% Export new solution from estimation
[sol_out,~]  = MapSolution(Wp,sol_out,Inf,options); % Map solution to flowfields
[~,sol_out]  = Actuator(Wp,sol_out,options);        % Recalculate power after analysis update
end