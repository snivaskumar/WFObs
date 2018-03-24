function [Wp,sol_out,strucObs] = WFObs_o_dexkf(strucObs,Wp,sys_in,sol_in,options)

if sol_in.k == 1
    if strucObs.fusion_type == 0     % CI = 0,1; EI = 2; ICI = 3, IFAC = 4
        disp('Type of fusion: CI (Naive Method)')
    elseif strucObs.fusion_type == 1
        disp('Type of fusion: CI')
    elseif strucObs.fusion_type == 2
        disp('Type of fusion: EI')
    elseif strucObs.fusion_type == 3
        disp('Type of fusion: ICI')   
    else
        disp('Type of fusion: IFAC')
    end
    if strucObs.typeCZ == 1         % 1 if Z = Co-Variance, 2 if Z = Information
        disp('Type of filter: Conv. KF')
    else
        disp('Type of filter: Info. KF')
    end
    if strucObs.Subsys_length <= 5
        fprintf('Subsystem Length: %.2d\n',strucObs.Subsys_length*Wp.turbine.Drotor);
    else
        fprintf('Subsystem Length: %d\n',strucObs.Subsys_length);
    end
end

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
    tur = Wp.turbine.N;
    if strucObs.stateEst || strucObs.measFlow
        stateLocArray = zeros(strucObs.size_output,2);
        for iii = 1:strucObs.size_output
            [~,loci,~]           = WFObs_s_sensors_nr2grid(iii,Wp.mesh);
            stateLocArray(iii,:) = [loci.x, loci.y];
        end
    end
    turbLocArray = zeros(Wp.turbine.N,2);
    for iii = 1:Wp.turbine.N
        turbLocArray(iii,:) = [Wp.turbine.Crx(iii),Wp.turbine.Cry(iii)];
    end
    strucObs.state = stateLocArray;
end

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
Ck          = strucObs.Htt;
[rC cC]     = size(Ck);

% ExKF forecast update
soltemp     = sol_in;
soltemp.k   = soltemp.k - 1;
[solf,sysf]	= WFSim_timestepping( soltemp, sys_in, Wp, options );       % Forward propagation

% tic
if (sol_in.k == 1) 
%     || (rem(sol_in.k,50) == 0)
    clear Fk Bk
    Fk(sysf.pRCM,sysf.pRCM) = sysf.A(sysf.pRCM,sysf.pRCM)\sysf.Al(sysf.pRCM,sysf.pRCM); % Linearized A-matrix at time k
    Bk(sysf.pRCM,:)         = sysf.A(sysf.pRCM,sysf.pRCM)\sysf.Bl(sysf.pRCM,:);         % Linearized B-matrix at time k
    if ~options.exportPressures 
        Fk     = Fk(1:strucObs.size_output,1:strucObs.size_output);
        Bk     = Bk(1:strucObs.size_output,:);
    end
%     aa  = double(abs(Fk)>1e-4);
%     p   = symrcm(aa);
%      
%     A       = Fk(p,p);
%     Bk      = Bk(p,:);
%     Ck      = Ck(:,p);
%     state   = stateLocArray(p,:);
%     
%     strucObs.p      = p;
%     strucObs.Fk     = A;
%     strucObs.Bk     = Bk;
%     strucObs.Ck     = Ck;
%     strucObs.state  = state;
    
    strucObs.Fk     = Fk;
    strucObs.Bk     = Bk;
    strucObs.Ck     = Ck;
%     strucObs.state  = stateLocArray;
end
% toc
Fk      = strucObs.Fk;
Bk      = strucObs.Bk;
Ck      = strucObs.Ck;
state   = strucObs.state;
% p       = strucObs.p;

n = length(Fk);
p = [1:n];

[rB cB]     = size(Bk);
Dk          = zeros(rC,cB);

% Neglect pressure terms
if ~options.exportPressures 
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

% tic
RD              = Wp.turbine.Drotor;
Subsys_length   = strucObs.Subsys_length;
type            = strucObs.fusion_type;
typeCZ          = strucObs.typeCZ;
if (sol_in.k == 1) 
%     || (rem(sol_in.k,50) == 0)
    [x,d,p, F,D,G,H,Q,R,l,n,x_est,x_unest, P_unest] = subsystem_turbine(sol_in, p,Fk,Bk,Ck,QQ,RR, tur,state,turbLocArray, Subsys_length,RD, Sk1k1);
    strucObs.subsystem.x = x;           strucObs.subsystem.d = d;
    strucObs.subsystem.F = F;           strucObs.subsystem.D = D;
    strucObs.subsystem.G = G;           strucObs.subsystem.H = H;
    strucObs.subsystem.Q = Q;           strucObs.subsystem.R = R;
    strucObs.subsystem.l = l;           strucObs.subsystem.n = n;
    strucObs.subsystem.x_est = x_est;   strucObs.subsystem.x_unest = x_unest;
    strucObs.subsystem.P_unest = P_unest;
end
x       = strucObs.subsystem.x;         d       = strucObs.subsystem.d;
F       = strucObs.subsystem.F;         D       = strucObs.subsystem.D;
G       = strucObs.subsystem.G;         H       = strucObs.subsystem.H;
Q       = strucObs.subsystem.Q;         R       = strucObs.subsystem.R;
l       = strucObs.subsystem.l;         n       = strucObs.subsystem.n;
x_est   = strucObs.subsystem.x_est;     x_unest = strucObs.subsystem.x_unest;
P_unest = strucObs.subsystem.P_unest;

[xkk Pkk] = distributed_linear( x,d,p,l,n, F,D,G,H,Q,R, y, xkk1,xk1k1,Sk1k1, x_est,x_unest, P_unest, type,typeCZ );
% toc

sol_out.x   = xkk;
strucObs.Pk = Pkk; 
% Export new solution from estimation
[sol_out,~]  = MapSolution(Wp,sol_out,Inf,options); % Map solution to flowfields
[~,sol_out]  = Actuator(Wp,sol_out,options);        % Recalculate power after analysis update
end