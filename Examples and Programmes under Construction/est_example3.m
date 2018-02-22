clear all
close all
clc

load('model');
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
% A       = sys2.A;
B       = sys2.B;
% C       = sys2.C;
C = [1,1,1,0,0,0;
     1,0,0,1,0,0;
     0,0,0,1,1,0;
     0,1,1,0,0,1];
 
% C = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]

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

%%

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

l = N;

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

%% Distributed Kalman Filter

xt1t1   = x01;
Pt1t1   = eye(n,n);
% u1      = u1';
% ynoisy  = ynoisy';
xtt     = zeros(n,600);
xtt(:,1)= xt1t1;

hr = 4;

x       = cell(hr,1);
x{1}    = [1,2,3]';
x{2}    = [1,4]';
x{3}    = [4,5]';
x{4}    = [2,3,6]';

d       = cell(hr,1);
d{1}    = [4,5,6]';
d{2}    = [2,3,5,6]';
d{3}    = [1,2,3,6]';
d{4}    = [1,4,5]';

F = cell(hr,1);
D = cell(hr,1);
G = cell(hr,1);
H = cell(hr,1);
Q = cell(hr,1);
R = cell(hr,1);
Ptt1 = cell(hr,1);

for k = 1:N
%         xtt1        = A*xt1t1 + B1*u1(:,k);
%         Ptmp = zeros(n,n);
    for i = 1:hr
        F{i} = A(x{i},x{i});
        D{i} = A(x{i},d{i});
        G{i} = B1( x{i},: );
        H{i} = C( i,: );
        Q{i} = QQ( x{i},x{i} );
        R{i} = RR( i,i );
        nl   = size(F{i});
        
        Slff{i} = Pt1t1( x{i},x{i} );            % Sff(k-1/k-1)_l
        Slfd{i} = Pt1t1( x{i},d{i} );            % Sfd(k-1/k-1)_l
        Sldd{i} = Pt1t1( d{i},d{i} );
        
        xtt1{i} = F{i}*xt1t1(x{i}) + D{i}*xt1t1(d{i}) + G{i}*u1(:,k);
        Ptt1{i} = F{i}*Slff{i}*F{i}' ...
                    + F{i}*Slfd{i}*D{i}' ...
                    + (F{i}*Slfd{i}*D{i}')' ...
                    + D{i}*Sldd{i}*D{i}' + Q{i};
        
%         Ptmp(x{i},x{i}) = Ptt1{i};
        Zlkk1{i} = inv(Ptt1{i});
    end
%     Zlkk1 = bandinv(Ptt1,1,x,hr,n,0);
    for i = 1:hr
        zlkk1{i}        = Zlkk1{i}*xtt1{i};
        ik              = H{i}(:,x{i})'*inv(R{i})*ynoisy(i,k);
        Ik              = H{i}(:,x{i})'*inv(R{i})*H{i}(:,x{i});
        ztt{i}          = zlkk1{i} + ik;
        Zt1t1{i}        = Zlkk1{i} + Ik;
    end
    z = cell(hr,1);
    Z = cell(hr,1);
    for i = 1:hr
        z{i} = ztt{i};
        Z{i} = Zt1t1{i};
    end
%     xf = x{1};
%     zf = z{1};
%     Zf = Z{1};
%     for i = 1:hr-1
%         xtmp = union(xf,x{i+1});
%         l = length(xtmp);
%         x_a = zeros(l,1);
%         x_b = zeros(l,1);
%         X_a = zeros(l,l);
%         X_b = zeros(l,l);
% 
%         x_a(xf) = zf;
%         x_b(x{i+1}) = z{i+1};
%         X_a(xf,xf) = Zf;
%         X_b(x{i+1},x{i+1}) = Z{i+1};
%         
% %         [Si,Di] = eig(X_a);
% %         [Sj,Dj] = eig(pinv(Di^0.5)*pinv(Si)*X_b*Si*pinv(Di^0.5));
% % 
% %         Dij = zeros(size(Dj));
% %         for ii = 1: length(Dij)
% %             Dij(ii,ii) = min(1,Dj(ii,ii));
% %         end
% % 
% %         Xij = Si*(Di^0.5)*Sj*Dij*pinv(Sj)*(Di^0.5)*pinv(Si);
% 
% %         w1 = X_a*pinv(X_a + X_b);
% %         w2 = X_b*pinv(X_a + X_b);
% %         xij = w1*x_b + w2*x_a;
% %                 
% %         zf  = x_a + x_b - xij;
% %         Zf  = X_a + X_b - Xij;
% 
%         ZA      = X_a;
%         ZB      = X_b;
%         ff      = @(w) trace(-(ZA + ZB - ZA*pinv(w*ZB + (1-w)*ZA)*ZB)); % min -f = max f
%         omega   = fminbnd(ff,0,1,optimset('Display','off'));
%         Xij     = ZA*pinv(omega*ZB + (1-omega)*ZA)*ZB;
%         Zf      = X_a + X_b - Xij;
%         
%         K = ( X_a - (omega)*Xij )*pinv(Zf);
%         L = ( X_b - (1 - omega)*Xij )*pinv(Zf);
%         zf = K*x_a + L*x_b;
%         xf  = xtmp;
%     end
%     Pt1t1   = pinv(Zf);
%     xt1t1   = Pt1t1*zf;
    
    [xt1t1,Pt1t1] = fuze(z,Z,x,hr,n);
    
    xtt(:,k)    = xt1t1;
end

figure
subplot(221)
set(gca,'FontSize',18)
plot(t,real(xtt(1,:)'),'b','LineWidth',2)
hold
plot(t,xnf(:,1),'r')
legend('Estimated','True')
title('Distributed Information KF - x1')
hold off

subplot(222)
plot(t,real(xtt(2,:)'),'b','LineWidth',2)
hold
plot(t,xnf(:,2),'r')
legend('Estimated','True')
title('Distributed Information KF - x2')
hold off

subplot(223)
set(gca,'FontSize',18)
plot(t,real(xtt(3,:)'),'b','LineWidth',2)
hold
plot(t,xnf(:,3),'r')
legend('Estimated','True')
title('Distributed Information KF - x1')
hold off

subplot(224)
plot(t,real(xtt(4,:)'),'b','LineWidth',2)
hold
plot(t,xnf(:,4),'r')
legend('Estimated','True')
title('Distributed Information KF - x2')
hold off
   