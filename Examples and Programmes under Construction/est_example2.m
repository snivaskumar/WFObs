clear all
close all
clc

load('model.mat');
L = 2;
for i = 1:6
    for j = 1:6
        if abs(i-j) <= L
            A(i,j) = sys2.A(i,j);
        else
            A(i,j) = 0;
        end
    end
end
B       = sys2.B;
C = [1,1,1,0,0,0;
     1,0,0,1,0,0;
     0,0,0,1,1,0;
     0,1,1,0,0,1];

% C = [1,1,1,1,0,0;
%      0,1,1,1,1,0;
%      0,0,1,1,1,1;
%      1,1,0,0,1,1];

D       = sys2.D;
n       = 6;
ts      = 0.01;

tmp     = rand(n,n);
Bw      = 0.05*tmp/norm(tmp);
tmp     = rand(4,n);
Dw      = 0.05*tmp/norm(tmp);

u       = [ones(1/ts,1); 0.2*ones(1/ts,1)];
u       = [u;u;u];
u       = [u,u];
N       = length(u);
t       = [0:ts:(N - 1)*ts]';

% noise covariance matrix 
pncov   = 0.01;     % Co-Variance of the Process Noise
mncov   = 0.01;     % Co-Variance of the Measurement Noise
Rprim   = mncov*eye(4,4);
Q1      = 0.01*eye(n,n);
R       = Dw*Q1*Dw' + Rprim;
S       = Bw*Q1*Dw';
Q       = Bw*Q1*Bw';
vprim   = randn(600,4)*((R-Dw*Q*Dw')^0.5);
w       = sqrt(pncov)*randn(600,6);
u1      = [u w vprim];

B1      = [B Bw zeros(6,4)];
D1      = [D Dw eye(4,4)];
T       = ss(A,B,C,D,ts);
T1      = ss(A,B1,C,D1,ts);


x01=[1.5 1.5 1.5 1.5 1.5 1.5]';

[ynoisy,t1,xnoisy]      = lsim(T1,u1,t,x01);
[ynf,t,xnf]             = lsim(T,u,t,x01);

QQ = Q;
RR = R;

%% Centralized KF

xk1k1   = x01;
Pk1k1   = eye(n,n);
u1      = u1';
ynoisy  = ynoisy';
xkk     = zeros(n,600);
xkk(:,1)= xk1k1;
for k = 1:N
    xkk1        = A*xk1k1 + B1*u1(:,k);
    Pkk1        = A*Pk1k1*A' + Q;
    
    dy          = ynoisy(:,k) - C*xkk1;
    Pyy         = R + C*Pkk1*C';
    Pxy         = Pkk1*C';
    K           = Pxy*inv(Pyy);
    xk1k1       = xkk1 + K*dy;
    Pk1k1       = ( eye(n,n) - K*C )*Pkk1;
    
    xkk(:,k)    = xk1k1;
end

l = 100;

figure
subplot(221)
set(gca,'FontSize',18)
plot(t(1:l),xkk(1,1:l)','b','LineWidth',2)
hold
plot(t(1:l),xnf(1:l,1),'r')
legend('Estimated','True')
title('Conventional KF - x1')
hold off

subplot(222)
plot(t(1:l),xkk(2,1:l)','b','LineWidth',2)
hold
plot(t(1:l),xnf(1:l,2),'r')
legend('Estimated','True')
title('Conventional KF - x2')

subplot(223)
plot(t(1:l),xkk(3,1:l)','b','LineWidth',2)
hold
plot(t(1:l),xnf((1:l),3),'r')
legend('Estimated','True')
title('Conventional KF - x3')

subplot(224)
plot(t(1:l),xkk(4,1:l)','b','LineWidth',2)
hold
plot(t(1:l),xnf(1:l,4),'r')
legend('Estimated','True')
title('Conventional KF - x4')
hold off

%% Information Filter

Pk1k1   = eye(n,n);
zk1k1   = x01;
zkk     = zeros(n,600);
zkk(:,1)= zk1k1;
for k = 1:N
    xkk1        = A*zk1k1 + B1*u1(:,k);
%     Zk1k1       = inv(Pk1k1);
%     Zkk1        = pinv(Q) - pinv(Q)*A*pinv(Zk1k1 + A'*pinv(Q)*A)*A'*pinv(Q);
    Zkk1        = inv( A*Pk1k1*A' + Q );
    zkk1        = Zkk1*xkk1;
    
    ik          = C'*inv(R)*ynoisy(:,k);
    Ik          = C'*inv(R)*C;
    zk1k1       = zkk1 + ik;
    Zk1k1       = Zkk1 + Ik;
    
    Pk1k1       = inv( Zk1k1 );
    zk1k1       = Pk1k1*zk1k1;
    zkk(:,k)    = zk1k1;
end

figure
subplot(221)
set(gca,'FontSize',18)
plot(t(1:l),zkk(1,1:l)','b','LineWidth',2)
hold
plot(t(1:l),xnf(1:l,1),'r')
legend('Estimated','True')
title('Information KF - x1')
hold off

subplot(222)
plot(t(1:l),zkk(2,1:l)','b','LineWidth',2)
hold
plot(t(1:l),xnf(1:l,2),'r')
legend('Estimated','True')
title('Information KF - x2')

subplot(223)
plot(t(1:l),zkk(3,1:l)','b','LineWidth',2)
hold
plot(t(1:l),xnf((1:l),3),'r')
legend('Estimated','True')
title('Information KF - x3')

subplot(224)
plot(t(1:l),zkk(4,1:l)','b','LineWidth',2)
hold
plot(t(1:l),xnf(1:l,4),'r')
legend('Estimated','True')
title('Information KF - x4')
hold off

%% DIstributed Information Filter

hr = 4;

% x       = cell(hr,1);
% x{1}    = [1,2,3]';
% x{2}    = [1,4]';
% x{3}    = [4,5]';
% x{4}    = [2,3,6]';
% 
% d       = cell(hr,1);
% d{1}    = [4,5,6]';
% d{2}    = [2,3,5,6]';
% d{3}    = [1,2,3,6]';
% d{4}    = [1,4,5]';

zt      = zeros(n,600);
Pt1t1   = eye(n,n);
zt1t1   = x01;

F = cell(hr,1);
D = cell(hr,1);
G = cell(hr,1);
H = cell(hr,1);
Q = cell(hr,1);
R = cell(hr,1);

xtt1 = cell(hr,1);
ztt1 = cell(hr,1);
zlkk1 = cell(hr,1);
Zlkk1 = cell(hr,1);
Ptmp = cell(hr,1);

type = 3;       % CI = 1; EI = 2; ICI = 3
tt = 1;
l = tt;
for k = 1:tt
    xkk1 = A*zt1t1 + B1*u1(:,k);
    [zt1t1 Pt1t1]   = subsystem( A,B1,C,D, ynoisy(:,k), xkk1,zt1t1,Pt1t1, QQ,RR, type );
    zt(:,k)         = zt1t1;
end
zt1t1
zkk(:,tt)

figure
subplot(221)
set(gca,'FontSize',18)
plot(t(1:l),real(zt(1,1:l)'),'b','LineWidth',2)
hold
plot(t(1:l),xnf(1:l,1),'r')
legend('Estimated','True')
title('Distributed Information KF - x1')
hold off

subplot(222)
plot(t(1:l),real(zt(2,1:l)'),'b','LineWidth',2)
hold
plot(t(1:l),xnf(1:l,2),'r')
legend('Estimated','True')
title('Distributed Information KF - x2')
hold off

subplot(223)
set(gca,'FontSize',18)
plot(t(1:l),real(zt(3,1:l)'),'b','LineWidth',2)
hold
plot(t(1:l),xnf(1:l,3),'r')
legend('Estimated','True')
title('Distributed Information KF - x1')
hold off

subplot(224)
plot(t(1:l),real(zt(4,1:l)'),'b','LineWidth',2)
hold
plot(t(1:l),xnf(1:l,4),'r')
legend('Estimated','True')
title('Distributed Information KF - x2')
hold off
