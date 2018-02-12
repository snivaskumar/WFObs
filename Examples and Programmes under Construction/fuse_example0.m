% clear all
% close all
% clc
% 
% hr = 4;
% x       = cell(hr,1);
% x{1}    = [1,2,3]';
% x{2}    = [1,4]';
% x{3}    = [4,5]';
% x{4}    = [2,3,6]';
% 
% zk1k1= [1.9278; -0.1936; 2.7214; 1.1245; 1.8360; 1.9433]
% load('/Users/Nivas_Kumar/Desktop/haha.mat','Zt1t1','ztt');
% Slkk = Zt1t1;
% xlkk = ztt;
% 
% i = 2;
% j = 1;
% a = xlkk{i};
% b = xlkk{j};
% P1 = Slkk{i};
% P2 = Slkk{j};
% 
% l1 = length(x{i});
% [cc1,ia1]=setdiff(x{i},x{j});
% Tmp1 = zeros(l1,l1);
% Tmp1(ia1,ia1) = eye(length(ia1),length(ia1));
% Tmp1 = Tmp1(ia1,[1:l1])
% 
% l2 = length(x{j});
% [ha,cc2,ia2] = intersect(x{i},x{j});
% Tmp2 = zeros(l2,l2);
% Tmp2(ia2,ia2) = eye(length(ia2),length(ia2));
% Tmp2 = Tmp2(ia2,[1:l2])
% 
% x_a = x{i}
% x_b = [Tmp2*x{j}; Tmp1*x{i}];
% 
% [x_b,ix] = sort(x_b);
% x_b
% 
% a       = a;
% btmp    = [Tmp2*b; Tmp1*a];
% b       = btmp(ix);
% % x_b
% % 
% lc1 = length(cc1);
% lc2 = length(cc2);
% 
% PP = zeros(l1,l1);
% PP(cc2,cc2) = Tmp2*P2*Tmp2';
% 
% [Sa,Da] = eig(P1);
% [Sb,Db] = eig(PP);
% 
% [~,I] = sort(diag(Db),'descend');
% Sb = Sb(I,I);
% Db = Db(I,I);
% 
% for ii = 1:length(x{i})
%     if Db(ii,ii) == 0
%         Gb(ii,ii) = 1;
%     else
%         Gb(ii,ii) = Db(ii,ii);
%     end
% end
% X_b = Sa*(Da^0.5)*Sb*Gb*pinv(Sb)*(Da^0.5)*pinv(Sa);
% X_a = P1;

%%
clear all
close all
clc

hr = 4;
x       = cell(hr,1);
x{1}    = [1,2,3]';
x{2}    = [1,4]';
x{3}    = [4,5]';
x{4}    = [2,3,6]';

zk1k1= [1.9278; -0.1936; 2.7214; 1.1245; 1.8360; 1.9433]
load('/Users/Nivas_Kumar/Desktop/haha.mat','Zt1t1','ztt');
Slkk = Zt1t1;
xlkk = ztt;

xf      = x{1};
zf      = xlkk{1};
Zf      = Slkk{1};
type    = 3;
for i = 1:hr-1
    x2  = x{i + 1}; 
    b   = xlkk{i + 1};
    P2  = Slkk{i + 1};

x_tmp                           = union( xf,x2 );
[cc1,ia1]                       = setdiff( x_tmp,xf );
T_tmp1                          = zeros(max(x_tmp),max(x_tmp));
T_tmp1(x_tmp,x_tmp)             = eye(length(x_tmp),length(x_tmp));
T_tmp1(cc1,cc1)                 = zeros(length(cc1),length(cc1));
T_tmp1                          = T_tmp1(x_tmp,xf);
T_tmp1                          = sparse(T_tmp1);
x_a                             = T_tmp1*zf;

[cc2,ia2]                       = setdiff( x_tmp,x2 );
T_tmp2                          = zeros(max(x_tmp),max(x_tmp));
T_tmp2(x_tmp,x_tmp)             = eye(length(x_tmp),length(x_tmp));
T_tmp2(cc2,cc2)                 = zeros(length(cc2),length(cc2));  
T_tmp2                          = T_tmp2(x_tmp,x2);
T_tmp2                          = sparse(T_tmp2);
x_b                             = T_tmp2*b;

X_a = T_tmp1*Zf*T_tmp1';
X_b = T_tmp2*P2*T_tmp2';

%%%%%%
% Extra
x_a(ia1) = x_b(ia1);
x_b(ia2) = x_a(ia2);

X_a(ia1,ia1) = X_b(ia1,ia1);
X_b(ia2,ia2) = X_a(ia2,ia2);
%%%%%%

if type == 0
    Zf      = 0.5.*X_a + 0.5.*X_b;
    zf      = 0.5.*x_a + 0.5.*x_b;
    0;
elseif type == 1
    % CI
    ZA      = X_a;
    ZB      = X_b;
    f       = @(w) trace( pinv(w*ZA + (1-w)*ZB) ); % arg (min -f) = arg (max f)
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
xf  = x_tmp;
end
P = pinv(Zf);
xk = P*zf

