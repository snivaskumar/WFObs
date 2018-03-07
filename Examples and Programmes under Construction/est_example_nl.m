clear all
close all
clc

n       = 1;    % Number of states
l       = 1;    % Number of Output
ts      = 1;    % Sampling Period
N       = 600;   % Estimation length
Ntmp    = N*ts;

tmp     = rand(n,n);
Bw      = 0.05*tmp/norm(tmp);
tmp     = rand(l,n);
Dw      = 0.05*tmp/norm(tmp);

utmp    = [ones(1/ts,1); 0.2*ones(1/ts,1)];
size(utmp);
u       = [];
for i = 1:(Ntmp/2)
    u       = [u;utmp];
end
utmp    = u;
size(utmp);
u       = [];
for i = 1:l
    u       = [u,utmp];
end
size(u);
t       = [0:ts:(N - 1)*ts]';

% noise covariance matrix 
pncov   = 0.01;     % Co-Variance of the Process Noise
mncov   = 0.01;     % Co-Variance of the Measurement Noise
Rprim   = mncov*eye(l,l);
Q1      = 0.01*eye(n,n);
R       = Dw*Q1*Dw' + Rprim;
S       = Bw*Q1*Dw';
Q       = Bw*Q1*Bw';
vprim   = randn(N,l)*((R-Dw*Q*Dw')^0.5);
w       = sqrt(pncov)*randn(N,n);
u1      = [u w vprim];

omega   = 4e-2;
a1      = 0.5; 
a2      = 0.2;
a3      = 0.5;
x01     = 1.5;

f = @(k,x,u,w,v) (1 + sin(pi*omega*k) + a1*x + w);
g1 = @(k,x,u,w,v) (a2.*(x.^2) + v);
g2 = @(k,x,u,w,v) (a3.*x - 2 + v);
x = zeros(1,600);
x(:,1) = x01;
for k = 1:N
    x(:,k+1) = f(k,x(:,k),0,w(k,:)',vprim(k,:)');
%     if k <= (N/2)
        y(:,k) = g1(k,x(:,k),0,w(k,:)',vprim(k,:)');
%     else
%         y(:,k) = g2(x(:,k),vprim(k,:)');
%     end
end
% plot(y')

%% ExKF

xk1k1   = x01;
Pk1k1   = eye(n,n);
u1      = u1';
xkk     = zeros(n,N);
xkk(:,1)= xk1k1;
for k = 1:N
%     xkk1        = A*xk1k1 + B1*u1(:,k);
%     Pkk1        = A*Pk1k1*A' + Q;
    
    A           = dfdx(f,k,xk1k1,0,w(k,:)',vprim(k,:)',n);
    C           = dfdx(g1,k,xk1k1,0,w(k,:)',vprim(k,:)',n);

    xkk1        = f(k,xk1k1,0,w(k,:)',vprim(k,:)');
    Pkk1        = A*Pk1k1*A' + Q;
    
%     dy          = ynoisy(:,k) - C*xkk1;
    dy          = y(:,k) - g1(k,xkk1,0,w(k,:)',vprim(k,:)');
    Pyy         = R + C*Pkk1*C';
    Pxy         = Pkk1*C';
    K           = Pxy*inv(Pyy);
    xk1k1       = xkk1 + K*dy;
    Pk1k1       = ( eye(n,n) - K*C )*Pkk1;
    
    xkk(:,k)    = xk1k1;
end
% figure,plot(xkk)

figure
set(gca,'FontSize',18)
plot(t,xkk','b','LineWidth',2)
hold
plot(t,y','r','LineWidth',2)
legend('Estimated','True')
title('ExKF')
hold off

%% UKF

%w = (1/5)*ones(1,5);
alpha   = .83;
beta    = 2;
k       = 0.0;
% alpha   = 1;
% beta    = 0;
% k       = .8;

lw      = length(w(1,:));
lv      = length(vprim(1,:));
L       = (n+lw+lv);

lambda  = (alpha^2)*(L + k) - L;
for i = 1:(2*L+1)
    if i == 1
        Wm(i)    = lambda/(L + lambda);
        Wc(i)    = Wm + (1 - alpha^2 + beta);
    else 
        Wm(i)    = 0.5/(L + lambda);
        Wc(i)    = Wm(i);
    end
end
Wm = Wm';
temp = eye(2*L+1) - repmat(Wm,1,2*L+1);
W    = temp*diag(Wc)*temp';

xt1t1   = x01;
Pt1t1   = eye(n,n);
u1      = u1';
xtt     = zeros(n,N);
xtt(:,1)= xt1t1;
for k = 1:N
    x_aug           = [xt1t1; w(k,:)'; vprim(k,:)'];
    P_aug           = [Pt1t1, zeros(n,lw+lv); zeros(n,lw+lv)' [R S; S Q]];
    
    ux              = repmat(x_aug,1,2*L+1); 
    Uscented_devs   = sqrt(L+lambda)*(sqrt(P_aug));
    ux(:,2:L+1)     = ux(:,2:L+1) + Uscented_devs;
    ux(:,L+2:end)   = ux(:,L+2:end) - Uscented_devs;
    
    xtt1        = f(k,ux(1:n,:),0,ux(n+1:n+lw,:),ux(n+lw+1:L,:));
    yt          = g1(k,ux(1:n,:),0,ux(n+1:n+lw,:),ux(n+lw+1:L,:));
    
    xtt1_aug    = [xtt1; ux(n+1:end,:)];
    yt_aug      = [yt; ux(n+1:end,:)];
    
    xmean       = sum(repmat(Wm',n,1) .* xtt1_aug, 2);
    ymean       = sum(repmat(Wm',l,1) .* yt_aug, 2);
    Aenft       = xtt1_aug - repmat(xmean,1,2*L+1);
    Yenft       = yt_aug - repmat(ymean,1,2*L+1);
    
    Ptt1        = Aenft*W*Aenft';
    
    dy          = [y(:,k);w(k,:)'; vprim(k,:)'] - ymean;
    Pyy         = Yenft*W*Yenft';
    Pxy         = Aenft*W*Yenft';
    K           = Pxy*pinv(Pyy);
    xt1t1_aug   = xmean + K*dy;
    Pt1t1_aug   = Ptt1 - K*Pyy*K';
    
    xt1t1       = xt1t1_aug(1:n);
    Pt1t1       = Pt1t1_aug(1:n,1:n);
    xtt(:,k)    = xt1t1;
end

figure
set(gca,'FontSize',18)
plot(t,xtt','b','LineWidth',2)
hold
plot(t,y','r','LineWidth',2)
legend('Estimated','True')
title('UKF')
hold off



