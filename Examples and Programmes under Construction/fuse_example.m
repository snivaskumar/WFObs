clear all
close all
clc

% a = rand(5,1);
% b = a;
% P = rand(5,5);
% P = (P*P')/norm(P*P');

a = [0.4480; 0.5611; 0.0698; 0.4981; 0.7489];
P = [0.3836    0.2449    0.1045    0.3815    0.1293;
     0.2449    0.1616    0.0665    0.2476    0.0809;
     0.1045    0.0665    0.0324    0.1066    0.0317;
     0.3815    0.2476    0.1066    0.3949    0.1252;
     0.1293    0.0809    0.0317    0.1252    0.0556];

Z = inv(P);
z = Z*a;


x = cell(3,1);
x{1} = [1,2,3]';
x{2} = [2,3,4]';
x{3} = [4,5]';

a1 = a(x{1});a2 = a(x{2});a3 = a(x{3});
z1 = z(x{1});z2 = z(x{2});z3 = z(x{3});
P1 = P(x{1},x{1});P2 = P(x{2},x{2});P3 = P(x{3},x{3});
Z1 = Z(x{1},x{1});Z2 = Z(x{2},x{2});Z3 = Z(x{3},x{3});

[c,C,cc] = fuze2(a1,a2,P1,P2,x{1},x{2},1);
[zc,ZC,zcc] = haha(z1,z2,Z1,Z2,x{1},x{2},3);
zc,z(zcc),ZC,Z(zcc,zcc)

function [zf, Zf, xf] = haha(z1,z2,Z1,Z2,x1,x2,type);

x_tmp                           = union( x1,x2 );
[cc{1},ia{1}]                   = setdiff( x_tmp,x1 );
T_tmp{1}                        = zeros(max(x_tmp),max(x_tmp));
T_tmp{1}(x_tmp,x_tmp)           = eye(length(x_tmp),length(x_tmp));
T_tmp{1}(cc{1},cc{1})           = zeros(length(cc{1}),length(cc{1}));
T_tmp{1}                        = T_tmp{1}(x_tmp,x1);
T_tmp{1}                        = sparse(T_tmp{1});
x_a                             = T_tmp{1}*z1;
x1 = T_tmp{1}*x1;

[cc{2},ia{2}]                   = setdiff( x_tmp,x2 );
T_tmp{2}                        = zeros(max(x_tmp),max(x_tmp));
T_tmp{2}(x_tmp,x_tmp)           = eye(length(x_tmp),length(x_tmp));
T_tmp{2}(cc{2},cc{2})           = zeros(length(cc{2}),length(cc{2}));  
T_tmp{2}                        = T_tmp{2}(x_tmp,x2);
T_tmp{2}                        = sparse(T_tmp{2});
x_b                             = T_tmp{2}*z2;
x2 = T_tmp{2}*x2;

X_a = T_tmp{1}*Z1*T_tmp{1}';
X_b = T_tmp{2}*Z2*T_tmp{2}';

%%%%%%
% Extra
x_a(ia{1}) = x_b(ia{1});
x_b(ia{2}) = x_a(ia{2});

x1(ia{1}) = x2(ia{1});
x2(ia{2}) = x1(ia{2});
x1;
x2;

X_a(ia{1},ia{1}) = X_b(ia{1},ia{1});
X_b(ia{2},ia{2}) = X_a(ia{2},ia{2});
%%%%%%

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
    ff      = @(w) trace(-(ZA + ZB - ZA*inv(w*ZB + (1-w)*ZA)*ZB)); % min -f = max f
    omega   = fminbnd(ff,0,1,optimset('Display','off'));
    Xij     = ZA*inv(omega*ZB + (1-omega)*ZA)*ZB;
    Zf      = X_a + X_b - Xij;

    K = ( X_a - (omega)*Xij )*pinv(Zf);
    L = ( X_b - (1 - omega)*Xij )*pinv(Zf);
    zf = K*x_a + L*x_b;
    3;
end
xf  = x_tmp;
end