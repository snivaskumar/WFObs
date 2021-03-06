function [zf,Zf,x_tmp] = haha(a,b,P1,P2,x1,x2,type);

a1 = a;
b2 = b;

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

[xx,x_c] = setdiff(x2,ha);
l = length(x_c);
c = b(x_c);
C = P2(x_c,x_c);

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
% [~, I] = sort(diag(Gb),'descend');
% Gb = Gb(I,I);
% Sb = Sb(I,I);

X_b = Sa*(Da^0.5)*Sb*Gb*pinv(Sb)*(Da^0.5)*pinv(Sa);
X_a = P1;

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

[qq,ww,ee] = intersect( x_tmp,ha );
DA = ZA(ww,ww);
DB = ZB(ww,ww);

[hatmp,ctmp,iatmp] = intersect( x_tmp,x2 );
x_a(ia1)    = x2(x_c);
x_b         = x_a;
a(ia1)      = b2(x_c);
b(ia1)      = b2(x_c);
% ZA(ia1,ia1) = P2(x_c,x_c);
% ZB(ia1,ia1) = P2(x_c,x_c);
ZA(ctmp,ctmp) = P2;
ZB(ctmp,ctmp) = P2;

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
% end