% Using observer and KF on the double tank
%
clear all
close all
clc
Ts=0.1; % sample time
N=600;  % the length of the simulation
n = 2;

% system matrixes
% A=[0.9512  0;
%    0.0476  0.9512];
A=[0.98  0;
   0.0076  0.9912];
B=[ 0.0975;
    0.0024];
Bw=[0.0975  0;
    0.0024  0.0975];
% C = [0 1];
C = [1 0];
rank(obsv(A,C))
Cd = [0 1; 1 0];
Dw=[0 0.5];
D=0;
T=ss(A,B,C,D,Ts);

% noise covariance matrix 
Rprim=0.01;
R1=0.0125;
Q1=[ 0.01    0;
    0       0.01];
S1=[0 0.005]';
R=Dw*Q1*Dw' + Rprim;
S=Bw*Q1*Dw';
Q=Bw*Q1*Bw';

% input signal
u=[ ones(10/0.1,1); 0.2*ones(10/0.1,1)];    
u=[u;u;u];
t=[0:0.1:599*0.1]';

% the noise signals
vprim=sqrt(R-Dw*Q*Dw')*randn(600,1);
w       =sqrt(0.01)*randn(600,2);
% wsrc=(Bw*w')';
% vsrc=(Dw*w'+vprim')';

% define a new system with input u1=[u w v]^T;
B1=[B Bw zeros(2,1)];
D1=[0 Dw 1];
Dd1=[0 Dw 1; 0 Dw 1];
T1=ss(A,B1,C,D1,Ts);
Td1=ss(A,B1,Cd,Dd1,Ts);
Tnf =ss(A,B,C,0,Ts);

% inital state
x01=[1.5 1.5]';

% the input
u1=[u w vprim];

% simulation
[ynoisy,t1,xnoisy]      = lsim(T1,u1,t,x01);
[ydnoisy,td1,xdnoisy]   = lsim(Td1,u1,t,x01);
% simulation noise free
[ynf,t1,xnf]=lsim(Tnf,u,t,x01);

f = @(x,u,w) (A*x + B*u + Bw*w);
g = @(x,u,w,v) (C*x + D*u + Dw*w + v);

% x = zeros(2,600);
% x(:,1) = x01;
% for k = 1:N
%     x(:,k+1) = f(x(:,k),u(k),w(k,:)');
%     y(:,k)   = g(x(:,k),u(k),w(k,:)',vprim(k,:)');
% end

%% Kalman Filter 

xk1k1   = x01;
Pk1k1   = eye(n,n);
u1      = u1';
ynoisy  = ynoisy';
xkk     = zeros(2,600);
xkk(:,1)= xk1k1;
% tmp     = cell(N,1); 
tmp1     = []; 
tmp2     = []; 
for k = 1:N
    xkk1        = A*xk1k1 + B1*u1(:,k);
    Pkk1        = A*Pk1k1*A' + Q;
    
    dy          = ynoisy(:,k) - C*xkk1;
    
    dyy(k)      = dy;
    
    Pyy         = R + C*Pkk1*C';
    Pxy         = Pkk1*C';
    K           = Pxy*inv(Pyy);
    xk1k1       = xkk1 + K*dy;
    Pk1k1       = ( eye(n,n) - K*C )*Pkk1;
    
    xkk(:,k)    = xk1k1;
    
    stmp1(k)      = max(abs(eig(Pkk1)));
    stmp2(k)      = max(abs(eig(Pk1k1)));
    
    tmp1(k)      = max(diag(Pkk1));
    tmp2(k)      = max(diag(Pk1k1));
end

XCsrcf=xcorr(dyy);
figure, plot(-floor(length(XCsrcf)/2):floor(length(XCsrcf)/2),XCsrcf,'LineWidth',2)
xlabel('Lag (\tau)')
axis([-100 100 -1 15])
hold on,
[Xps_srcf,om]=pwelch(dyy,128,[],[],1/Ts);
figure, semilogy(om,Xps_srcf,'LineWidth',2);
xlabel('\omega (Hz)')
set(gca,'FontSize',18)

figure
subplot(211)
set(gca,'FontSize',18)
plot(t,xkk(1,:)','b','LineWidth',2)
hold
plot(t,xnf(:,1),'r')
legend('Estimated','True')
title('Conventional KF - x1')
hold off

subplot(212)
plot(t,xkk(2,:)','b','LineWidth',2)
hold
plot(t,xnf(:,2),'r')
legend('Estimated','True')
title('Conventional KF - x2')
hold off

%% Information Filter

Pk1k1   = eye(n,n);
zk1k1   = x01;
zkk     = zeros(2,600);
zkk(:,1)= zk1k1;
for k = 1:N
    xkk1        = A*zk1k1 + B1*u1(:,k);
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
subplot(211)
set(gca,'FontSize',18)
plot(t,zkk(1,:)','b','LineWidth',2)
hold
plot(t,xnf(:,1),'r')
legend('Estimated','True')
title('Information KF - x1')
hold off

subplot(212)
plot(t,zkk(2,:)','b','LineWidth',2)
hold
plot(t,xnf(:,2),'r')
legend('Estimated','True')
title('Information KF - x2')
hold off

%% Distributed Information Filter

hr = 2;
ydnoisy  = ydnoisy';

x       = cell(hr,1);
x{1}    = [1,2]';
x{2}    = [1,2]';

Zt1t1   = cell(hr,1);
Pt1t1   = cell(hr,1);
ztt     = cell(hr,1);
zt      = zeros(2,600);
Pt1t1   = eye(n,n);
zt1t1   = x01;

for k = 1:N
    for i = 1:hr
        xtt1            = A*zt1t1 + B1*u1(:,k);
        Ztt1            = inv( A*Pt1t1*A' + Q );
        ztt1            = Ztt1*xtt1;

        ik              = Cd(1,:)'*inv(R)*ydnoisy(i,k);
        Ik              = Cd(1,:)'*inv(R)*Cd(1,:);
        ztt{i}          = ztt1 + ik;
        Zt1t1{i}        = Ztt1 + Ik;
    end
    Ptmp            = bandinv(Zt1t1,2,x,2,2,0);
%     [X, X_a, X_b, Xij, xxx, xab_size, T1, T2, ia1, ia2, cx] ...
%                 = fuse_cov(Ptmp{1},Ptmp{2},x{1},x{2});
%     ztt{1}          = Ptmp{1}*ztt{1};
%     ztt{2}          = Ptmp{2}*ztt{2};
%     [zt1t1, Pt1t1] = fuse_mean(ztt{1},ztt{2}, X, X_a, X_b, Xij, xxx, xab_size, T1, T2, ia1, ia2, cx);

    ztt{1}          = Ptmp{1}*ztt{1};
    ztt{2}          = Ptmp{2}*ztt{2};
    [zt1t1, Pt1t1]  = fuse2(ztt{1},ztt{2},Ptmp{1},Ptmp{2},x{1},x{2});
    zt(:,k)         = zt1t1;
end

% for k = 1:N
%     for i = 1:hr
%         xtt1            = A*zt1t1 + B1*u1(:,k);
%         Ztt1            = inv( A*Pt1t1*A' + Q );
%         ztt1            = Ztt1*xtt1;
% 
%         ik              = Cd(1,:)'*inv(R)*ydnoisy(i,k);
%         Ik              = Cd(1,:)'*inv(R)*Cd(1,:);
%         ztt{i}          = ztt1 + ik;
%         Zt1t1{i}        = Ztt1 + Ik;
%     end
%     
%     Ptmp            = bandinv(Zt1t1,2,x,2,2,0);
%     ztt{1}          = Ptmp{1}*ztt{1};
%     ztt{2}          = Ptmp{2}*ztt{2};
%     [zt1t1, Pt1t1]  = fuse2(ztt{1},ztt{2},Ptmp{1},Ptmp{2},x{1},x{2});
%     zt(:,k)         = zt1t1;
% end

figure
subplot(211)
set(gca,'FontSize',18)
plot(t,real(zt(1,:)'),'b','LineWidth',2)
hold
plot(t,xnf(:,1),'r')
legend('Estimated','True')
title('Distributed Information KF - x1')
hold off
subplot(212)
plot(t,real(zt(2,:)'),'b','LineWidth',2)
hold
plot(t,xnf(:,2),'r')
legend('Estimated','True')
title('Distributed Information KF - x2')
hold off



