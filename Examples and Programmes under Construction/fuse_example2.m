clear all
close all
clc

hr = 4;
x       = cell(hr,1);
x{1}    = [1,2,3]';
x{2}    = [1,4]';
x{3}    = [4,5]';
x{4}    = [2,3,6]';

zk1k1 = [1.8121; -0.1938; 2.7738; 1.2965; 1.7856; 1.7763];
load('/Users/Nivas_Kumar/Desktop/haha.mat','Zt1t1','ztt');
Slkk = Zt1t1;
xlkk = ztt;

i = 3;
j = 2;
a = xlkk{i};
b = xlkk{j};
P1 = Slkk{i};
P2 = Slkk{j};

l1 = length(x{i});
[cc1,ia1]=setdiff(x{i},x{j});
Tmp1 = zeros(l1,l1);
Tmp1(ia1,ia1) = eye(length(ia1),length(ia1));
Tmp1 = Tmp1(ia1,[1:l1]);

l2 = length(x{j});
[ha,cc2,ia2] = intersect(x{i},x{j});
Tmp2 = zeros(l2,l2);
Tmp2(ia2,ia2) = eye(length(ia2),length(ia2));
Tmp2 = Tmp2(ia2,[1:l2]);

[xx,x_c] = setdiff(x{j},ha);
l = length(x_c);
c = b(x_c);
C = P2(x_c,x_c);

x_a = x{i};
x_b = [Tmp1*x{i}; Tmp2*x{j}];

[x_b,ix] = sort(x_b);

a       = a;
btmp    = [Tmp1*a; Tmp2*b];
b       = btmp(ix);
 
lc1 = length(cc1);
lc2 = length(cc2);

PP = zeros(l1,l1);
PP(cc2,cc2) = Tmp2*P2*Tmp2';

[Sa,Da] = eig(P1);
[Sb,Db] = eig(PP);

for ii = 1:length(x{i})
    if Db(ii,ii) == 0
        Gb(ii,ii) = 1;
    else
        Gb(ii,ii) = Db(ii,ii);
    end
end
% [~, I] = sort(diag(Gb),'descend');
% Gb = Gb(I,I);
% Sb = Sb(I,I);

X_b = Sa*(Da^0.5)*Sb*Gb*pinv(Sb)*(Da^0.5)*pinv(Sa)
X_a = P1

la = length(a);
lb = length(b);

xa = [x_a; zeros(l,1)];
xb = [x_b; xx];
[xb,ix] = sort(xb);

a = [a; zeros(l,1)];
b = [b; c];
X_a = [X_a, zeros(la,l); zeros(l,la), zeros(l,l)];
X_b = [X_b, zeros(lb,l); zeros(l,lb), C];

a = a(ix)
b = b(ix)
X_a = X_a(ix,ix)
X_b = X_b(ix,ix)

ZA      = X_a;
ZB      = X_b;
x_a     = a;
x_b     = b;

ZA      = X_a;
ZB      = X_b;
f       = @(w) trace( -(w*ZA + (1-w)*ZB) ); % arg (min -f) = arg (max f)
omega   = fminbnd(f,0,1,optimset('Display','off'));
Zf      = omega*ZA + (1-omega)*ZB;
zf      = omega*x_a + (1-omega)*x_b;

% ff      = @(w) trace(-(ZA + ZB - ZA*pinv(w*ZB + (1-w)*ZA)*ZB)); % min -f = max f
% omega   = fminbnd(ff,0,1,optimset('Display','off'));
% Xij     = ZA*pinv(omega*ZB + (1-omega)*ZA)*ZB;
% Zf      = X_a + X_b - Xij;
% 
% K = ( X_a - (omega)*Xij )*pinv(Zf);
% L = ( X_b - (1 - omega)*Xij )*pinv(Zf);
% zf = K*x_a + L*x_b;

xf = inv(Zf)*zf

% [zff,Zff] = fuse(xlkk,Slkk,x,hr);
% inv(Zff)*zff