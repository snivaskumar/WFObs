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
% 
zk1k1= [1.9278; -0.1936; 2.7214; 1.1245; 1.8360; 1.9433]
load('haha.mat','Zt1t1','ztt');
Slkk = Zt1t1;
xlkk = ztt;

z = xlkk;
Z = Slkk;
z1 = z{1};z2 = z{2};z3 = z{3};z4 = z{4};
Z1 = Z{1};Z2 = Z{2};Z3 = Z{3};Z4 = Z{4};

% x{1} = [1,2,3]';
% x{2} = [2,3,4]';
% x{3} = [4,5]';
% a = [0.4480; 0.5611; 0.0698; 0.4981; 0.7489]
% P = [0.3836    0.2449    0.1045    0.3815    0.1293;
%      0.2449    0.1616    0.0665    0.2476    0.0809;
%      0.1045    0.0665    0.0324    0.1066    0.0317;
%      0.3815    0.2476    0.1066    0.3949    0.1252;
%      0.1293    0.0809    0.0317    0.1252    0.0556];
 
% a = rand(6,1)
% P = bandmatrix(6,1,1);
% P = rand(6,6);
% P = (P*P')/norm(P*P');
% Z = inv(P);
% z = Z*a;
% 
% a1 = a(x{1});a2 = a(x{2});a3 = a(x{3});a4 = a(x{4});
% z1 = z(x{1});z2 = z(x{2});z3 = z(x{3});z4 = z(x{4});
% P1 = P(x{1},x{1});P2 = P(x{2},x{2});P3 = P(x{3},x{3});P4 = P(x{4},x{4});
% Z1 = Z(x{1},x{1});Z2 = Z(x{2},x{2});Z3 = Z(x{3},x{3});Z4 = Z(x{4},x{4});
type = 3;
[zc,ZC,zcc] = haha(z1,z2,Z1,Z2,x{1},x{2},type);
[zc,ZC,zcc] = haha(zc,z3,ZC,Z3,zcc,x{3},type);
[zc,ZC,zcc] = haha(zc,z4,ZC,Z4,zcc,x{4},type);
real(pinv(ZC));
pinv(ZC)*zc
zc;
real(ZC);

[zc,ZC,zcc] = fuze2(z1,z2,Z1,Z2,x{1},x{2},type);
[zc,ZC,zcc] = fuze2(zc,z3,ZC,Z3,zcc,x{3},type);
[zc,ZC,zcc] = fuze2(zc,z4,ZC,Z4,zcc,x{4},type);
real(pinv(ZC));
pinv(ZC)*zc
zc;
real(ZC);

[x,C] = IFAC(z,Z,x,hr,n);
x

type = 3;
% i=2;j=3;
% [zf,Zf,xf] = haha(a{i},a{j},Z{i},Z{j},x{i},x{j},type);
% [zf,Zf,xf] = haha(zf,a{3},Zf,Z{3},xf,x{3},type);
% [zf,Zf,xf] = haha(zf,a{4},Zf,Z{4},xf,x{4},type);
% xf = pinv(Zf)*zf




function [a,b,X_a,X_b,x_a,x_b] = hello(a,b,P1,P2,x1,x2);

l1 = length(x1);
[cc1,ia1]=setdiff(x1,x2);
Tmp1 = zeros(l1,l1);
Tmp1(ia1,ia1) = eye(length(ia1),length(ia1));
Tmp1 = Tmp1(ia1,[1:l1]);

l2 = length(x2);
[ha,cc2,ia2] = intersect(x1,x2);
Tmp2 = zeros(l2,l2);
Tmp2(ia2,ia2) = eye(length(ia2),length(ia2));
Tmp2 = Tmp2(ia2,[1:l2]);

x_a = x1;
x_b = [Tmp1*x1; Tmp2*x2];

[x_b,ix] = sort(x_b);
x_b;
a       = a;
btmp    = [Tmp1*a; Tmp2*b];
b       = btmp(ix);

lc1 = length(cc1);
lc2 = length(cc2);

PP = zeros(l1,l1);
PP(cc2,cc2) = Tmp2*P2*Tmp2';

[Sa,Da] = eig(P1);
[Sb,Db] = eig(pinv(Da^0.5)*pinv(Sa)*PP*pinv(Sa')*pinv(Da^0.5));

for ii = 1:length(x1)
    if abs(Db(ii,ii)) <= 1e-4
        Gb(ii,ii) = 1;
    else
        Gb(ii,ii) = Db(ii,ii);
    end
end

X_b = Sa*(Da^0.5)*Sb*Gb*pinv(Sb)*(Da^0.5)*pinv(Sa);
X_a = P1;

end

function [zf,Zf,x_tmp] = haha(a,b,P1,P2,x1,x2,type);

a1 = a;
b2 = b;

[a,b,X_a,X_b,x_a,x_b] = hello(a1,b2,P1,P2,x1,x2);
[aa,bb,X_aa,X_bb,x_aa,x_bb] = hello(b2,a1,P2,P1,x2,x1);

x_tmp                           = union(x1,x2);
[cc1,ia1]                       = setdiff( x_tmp,x_a );
T_tmp1                          = zeros(max(x_tmp),max(x_tmp));
T_tmp1(x_tmp,x_tmp)             = eye(length(x_tmp),length(x_tmp));
T_tmp1(cc1,cc1)                 = zeros(length(cc1),length(cc1));
T_tmp1                          = T_tmp1(x_tmp,x_a);
T_tmp1                          = sparse(T_tmp1);
x_a                             = T_tmp1*x_a;

ZA      = X_a;
ZB      = X_b;

a   = T_tmp1*a;
b   = T_tmp1*b;
ZA  = T_tmp1*ZA*T_tmp1';
ZB  = T_tmp1*ZB*T_tmp1';

[ha,cc2,ia2] = intersect(x1,x2);
[xx,x_c] = setdiff(x2,ha);
[qq,ww,ee] = intersect( x_tmp,ha );
DA = ZA(ww,ww);
DB = ZB(ww,ww);

[hatmp,ctmp,iatmp] = intersect( x_tmp,x2 );
x_a(ia1)    = x2(x_c);
x_b         = x_a;
a(ia1)      = b2(x_c);
b(ia1)      = b2(x_c);
% ZA(ia1,ia1) = X_aa(x_c,x_c);
% ZB(ia1,ia1) = X_aa(x_c,x_c);
ZA(ctmp,ctmp) = X_aa;
ZB(ctmp,ctmp) = X_aa;

ZA(ww,ww) = DA;
ZB(ww,ww) = DB;

if type == 0
    Zf      = 0.5.*ZA + 0.5.*ZB;
    zf      = 0.5.*a + 0.5.*b;
    0;
elseif type == 1
    % CI
    f       = @(w) trace( -(w*ZA + (1-w)*ZB) ); % arg (min -f) = arg (max f)
    omega   = fminbnd(f,0,1,optimset('Display','off'));
    Zf      = omega*ZA + (1-omega)*ZB;
    zf      = omega*a + (1-omega)*b;
    1;
elseif type == 2
    % EI
    [Si,Di] = eig(ZA);
    [Sj,Dj] = eig(pinv(Di^0.5)*pinv(Si)*ZB*Si*pinv(Di^0.5));

    Dij = zeros(size(Dj));
    for ii = 1: length(Dij)
        Dij(ii,ii) = min(1,Dj(ii,ii));
    end
    Xij = Si*(Di^0.5)*Sj*Dij*pinv(Sj)*(Di^0.5)*pinv(Si);

%         w1 = (X_b - 0.5*Xij)*pinv(X_a + X_b - Xij);
%         w2 = (X_a - 0.5*Xij)*pinv(X_a + X_b - Xij);
        
    w1 = (ZB)*pinv(ZA + ZB);
    w2 = (ZA)*pinv(ZA + ZB);
    xij = w1*a + w2*b;

    zf  = a + b - xij;
    Zf  = ZA + ZA - Xij;
    2;
elseif type == 3    
    % ICI
    ff      = @(w) trace(-(ZA + ZB - ZA*pinv(w*ZB + (1-w)*ZA)*ZB)); % min -f = max f
    omega   = fminbnd(ff,0,1,optimset('Display','off'));
    Xij     = ZA*pinv(omega*ZB + (1-omega)*ZA)*ZB;
    Zf      = ZA + ZB - Xij;
    K = ( ZA - (omega)*Xij )*pinv(Zf);
    L = ( ZB - (1 - omega)*Xij )*pinv(Zf);
    zf = K*a + L*b;
    3;
end
end
