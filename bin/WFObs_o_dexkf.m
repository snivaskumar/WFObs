function [Wp,sol_out,strucObs] = WFObs_o_dexkf(strucObs,Wp,sys_in,sol_in,options)
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

if sol_in.k == 1
    xk1k1 = [vec(sol_in.u(3:end-1,2:end-1)'); vec(sol_in.v(2:end-1,3:end-1)')];
    if options.exportPressures == 1 % Optional: add pressure terms
        xk1k1 = [xk1k1; vec(sol_in.p(2:end-1,2:end-1)')];
        xk1k1 = xk1k1(1:end-2); % Correction for how pressure is formatted
    end
else
    xk1k1 = sol_in.x;
end
        
% Import measurement data variable
measuredData = sol_in.measuredData;

if sol_in.k == 1
    % Setup covariance and system output matrices
    if options.exportPressures
        strucObs.Pk    = sparse(eye(strucObs.size_state))*strucObs.P_0;
        strucObs.Htt   = sparse(eye(strucObs.size_state));
        strucObs.Htt   = strucObs.Htt(strucObs.obs_array,:);
        strucObs.Q_k   = strucObs.Q_k*eye(strucObs.size_state);
    else
        strucObs.Pk    = sparse(eye(strucObs.size_output))*strucObs.P_0;
        strucObs.Htt   = sparse(eye(strucObs.size_output));
        strucObs.Htt   = strucObs.Htt(strucObs.obs_array,:);
        strucObs.Q_k   = strucObs.Q_k*eye(strucObs.size_output);
    end;
end;

% ExKF forecast update
[solf,sysf]             = WFSim_timestepping( sol_in, sys_in, Wp, options ); % Forward propagation

if (sol_in.k == 1) || (rem(sol_in.k,5) == 0)
    clear Fk Bk
    Fk(sysf.pRCM,sysf.pRCM) = sysf.A(sysf.pRCM,sysf.pRCM)\sysf.Al(sysf.pRCM,sysf.pRCM); % Linearized A-matrix at time k
    Bk(sysf.pRCM,:)         = sysf.A(sysf.pRCM,sysf.pRCM)\sysf.Bl(sysf.pRCM,:);         % Linearized B-matrix at time k
    strucObs.Fk = Fk;
    strucObs.Bk = Bk;
end
Fk = strucObs.Fk;
Bk = strucObs.Bk;
[rB cB] = size(Bk);
Ck = strucObs.Htt;
[rC cC] = size(Ck);
Dk = zeros(rC,cB);

% Neglect pressure terms
if ~options.exportPressures 
    Fk     = Fk(1:strucObs.size_output,1:strucObs.size_output);
    Bk     = Bk(1:strucObs.size_output,:);
    solf.x = solf.x(1:strucObs.size_output);
end;

% ExKF analysis update
sol_out = sol_in;
xkk1    = solf.x;
y       = measuredData.sol(strucObs.obs_array);
Sk1k1   = strucObs.Pk;
QQ      = strucObs.Q_k;
lop     = length(strucObs.obs_array);
RR      = strucObs.R_k*eye(lop,lop);

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Centralised Information Filter
% tic
% Pf = Fk*strucObs.Pk*Fk' + strucObs.Q_k;  % Covariance matrix P for x(k) knowing y(k-1)
% Zkk1 = inv(Pf);
% zkk1 = Zkk1*xkk1;
% 
% ik          = Ck'*inv(RR)*y;
% Ik          = Ck'*inv(RR)*Ck;
% zkk         = zkk1 + ik;
% Zkk         = Zkk1 + Ik;
% 
% Pkk         = inv(Zkk);
% xkk         = Pkk*zkk;
% toc

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Distributed Information Filter
tic
type = 4;       % CI = 0,1; EI = 2; ICI = 3
[xkk Pkk]   = subsystem( Fk,Bk,Ck,Dk, y,xkk1,xk1k1,Sk1k1, QQ,RR, type );
toc

sol_out.x   = xkk;
strucObs.Pk = Pkk; 
% Export new solution from estimation
[sol_out,~]  = MapSolution(Wp,sol_out,Inf,options); % Map solution to flowfields
[~,sol_out]  = Actuator(Wp,sol_out,options);        % Recalculate power after analysis update
end