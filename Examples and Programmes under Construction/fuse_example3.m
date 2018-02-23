clear all
close all
clc

hr = 4;
n = 6;
x       = cell(hr,1);
x{1}    = [1,2,3]';
x{2}    = [1,4]';
x{3}    = [4,5]';
x{4}    = [2,3,6]';

zk1k1= [1.9278; -0.1936; 2.7214; 1.1245; 1.8360; 1.9433]
load('haha.mat','Zt1t1','ztt');
Slkk = Zt1t1;
xlkk = ztt;

z = ztt;
Z = Zt1t1;

xf = x{1};
zf = z{1};
Zf = Z{1};
type = 1
for i = 1:hr-1
    xtmp = union(xf,x{i+1});
    l = max(xtmp);
    x_a = zeros(l,1);
    x_b = zeros(l,1);
    X_a = zeros(l,l);
    X_b = zeros(l,l);
    
    x_a(xf) = zf(xf);
    x_b(x{i+1}) = z{i+1};
    X_a(xf,xf) = Zf(xf,xf);
    X_b(x{i+1},x{i+1}) = Z{i+1};
    
    x_a         = x_a(xtmp);
    x_b         = x_b(xtmp);
    X_a         = X_a(xtmp,xtmp);
    X_b         = X_b(xtmp,xtmp);
    
    [Si,Di] = eig(X_a);
    [Sj,Dj] = eig(pinv(sqrt(Di))*pinv(Si)*X_b*Si*pinv(sqrt(Di)));

    Dij = zeros(size(Dj));
    for ii = 1: length(Dij)
        Dij(ii,ii) = min(1,Dj(ii,ii));
    end
    
    if type == 0
        Zf      = X_a + X_b;
        zf      = x_a + x_b;
        0;
    elseif type == 1
        % CI
        ZA      = X_a;
        ZB      = X_b;
        f       = @(w) trace( -(w*ZA + (1-w)*ZB) ); % arg (min -f) = arg (max f)
        omega   = fminbnd(f,0,1,optimset('Display','off'));
        Zf      = omega*ZA + (1-omega)*ZB;
        zf      = omega*x_a + (1-omega)*x_b;
        1;
    elseif type == 2
        % EI
        [Si,Di] = eig(X_a);
        [Sj,Dj] = eig(pinv(Di^0.5)*pinv(Si)*X_b*Si*pinv(Di^0.5));

        Dij = zeros(size(Dj));
        for ii = 1: length(Dij)
            Dij(ii,ii) = min(1,Dj(ii,ii));
        end
        Xij = Si*(Di^0.5)*Sj*Dij*pinv(Sj)*(Di^0.5)*pinv(Si);

    %         w1 = (X_b - 0.5*Xij)*pinv(X_a + X_b - Xij);
    %         w2 = (X_a - 0.5*Xij)*pinv(X_a + X_b - Xij);

        w1 = (X_b)*pinv(X_a + X_b);
        w2 = (X_a)*pinv(X_a + X_b);
        xij = w1*x_a + w2*x_b;

        zf  = x_a + x_b - xij;
        Zf  = X_a + X_b - Xij;
        2;
    elseif type == 3    
        % ICI
        ZA      = X_a;
        ZB      = X_b;
        ff      = @(w) trace(-(ZA + ZB - ZA*pinv(w*ZB + (1-w)*ZA)*ZB)); % min -f = max f
        omega   = fminbnd(ff,0,1,optimset('Display','off'));
        Xij     = ZA*pinv(omega*ZB + (1-omega)*ZA)*ZB;
        Zf      = X_a + X_b - Xij;

        K = ( X_a - (omega)*Xij )*pinv(Zf);
        L = ( X_b - (1 - omega)*Xij )*pinv(Zf);
        zf = K*x_a + L*x_b;
        3;
    end
    xf  = xtmp;
end
P = pinv(Zf);
xk = P*zf

[zf, Zf, xf] = fuze2(ztt{1},ztt{2},Zt1t1{1},Zt1t1{2},x{1},x{2},type);
[zf, Zf, xf] = fuze2(zf,ztt{3},Zf,Zt1t1{3},xf,x{3},type);
[zf, Zf, xf] = fuze2(zf,ztt{4},Zf,Zt1t1{4},xf,x{4},type);
% real(Zf);
Zf = pinv(Zf);
% Zf;
xf = Zf*zf

[zf1, Zf1, xf1] = fuze2(ztt{1},ztt{2},Zt1t1{1},Zt1t1{2},x{1},x{2},type);
[zf2, Zf2, xf2] = fuze2(ztt{3},ztt{4},Zt1t1{3},Zt1t1{4},x{3},x{4},type);
[zff, Zff, xff] = fuze2(zf1,zf2,Zf1,Zf2,xf1,xf2,type);
pinv(Zff)*zff
real(Zff);
zff;

[xe,Ce] = fuze(z,Z,x,hr,n);
xe

% [zf, Zf, xf] = haha(ztt{1},ztt{2},Zt1t1{1},Zt1t1{2},x{1},x{2},type);
% [zf, Zf, xf] = haha(zf,ztt{3},Zf,Zt1t1{3},xf,x{3},type);
% [zf, Zf, xf] = haha(zf,ztt{4},Zf,Zt1t1{4},xf,x{4},type);
% real(Zf);
% zf
% Zf = pinv(Zf);
% % Zf;
% xf = Zf*zf
% 
% [zf1, Zf1, xf1] = haha(ztt{1},ztt{2},Zt1t1{1},Zt1t1{2},x{1},x{2},type);
% [zf2, Zf2, xf2] = haha(ztt{3},ztt{4},Zt1t1{3},Zt1t1{4},x{3},x{4},type);
% [zff, Zff, xff] = haha(zf1,zf2,Zf1,Zf2,xf1,xf2,type);
% pinv(Zff)*zff

%%

% clear all
% close all
% clc

% hr = 4;
% n = 6;
% x       = cell(hr,1);
% x{1}    = [1,2,3]';
% x{2}    = [1,4]';
% x{3}    = [4,5]';
% x{4}    = [2,3,6]';
% 
% zk1k1= [1.9278; -0.1936; 2.7214; 1.1245; 1.8360; 1.9433]

% load('haha.mat','Zt1t1','ztt');
% z = ztt;
% Z = Zt1t1;
% 
% % load('haha1.mat','Ptmp','ztt');
% % z = ztt;
% % Z = Ptmp;
% 
% H = [];
% k = 0;
% nn = 0;
% for i = 1:hr
%     l               = length(x{i});
%     nn              = nn + l;
%     H(k+1:nn,:)     = zeros(l,n);
%     H(k+1:nn,x{i})  = eye(l,l);
%     k               = k + l;
% end
% 
% C = [];
% nn = 0;
% k = 0;
% for i = 1:hr
%     l               = length(x{i});
%     nn              = nn + l;
%     CC              = inv(Z{i});     
% %     C(k+1:nn,k+1:nn)= (1/l).*CC;
%     C(k+1:nn,k+1:nn)= inv((1/l).*CC);
%     xx(k+1:nn,:)    = CC*z{i};
%     k               = k + l;
% end
% 
% % Ze = inv(H'*inv(C)*H);
% % ze = Ze*H'*inv(C)*xx
% 
% Ze = inv(H'*(C)*H);
% ze = Ze*H'*(C)*xx
% % Ce = inv(Ze);
% % xe = Ce*ze