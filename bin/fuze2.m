function [zf, Zf, xf] = fuze2(z1,z2,Z1,Z2,x1,x2,type);
% [zf, Zf, xf] = fuze2(z1,z2,Z1,Z2,x1,x2,type);
% type    1 for CI, 2 for EI, 3 for ICI

x_tmp                           = union( x1,x2 );
[cc{1},ia{1}]                   = setdiff( x_tmp,x1 );
T_tmp{1}                        = zeros(max(x_tmp),max(x_tmp));
T_tmp{1}(x_tmp,x_tmp)           = eye(length(x_tmp),length(x_tmp));
T_tmp{1}(cc{1},cc{1})           = zeros(length(cc{1}),length(cc{1}));
% T_tmp{i}                        = T_tmp{i}(x_tmp,setdiff(x_tmp,ia{i}));
T_tmp{1}                        = T_tmp{1}(x_tmp,x1);
T_tmp{1}                        = sparse(T_tmp{1});
x_a                             = T_tmp{1}*z1;
x1 = T_tmp{1}*x1;

% clear cc ia
[cc{2},ia{2}]                   = setdiff( x_tmp,x2 );
T_tmp{2}                        = zeros(max(x_tmp),max(x_tmp));
T_tmp{2}(x_tmp,x_tmp)           = eye(length(x_tmp),length(x_tmp));
T_tmp{2}(cc{2},cc{2})           = zeros(length(cc{2}),length(cc{2}));  
% T_tmp{i+1}                      = T_tmp{i+1}(x_tmp,setdiff(x_tmp,ia{i+1}));
T_tmp{2}                        = T_tmp{2}(x_tmp,x2);
T_tmp{2}                        = sparse(T_tmp{2});
x_b                             = T_tmp{2}*z2;
x2 = T_tmp{2}*x2;

X_a = T_tmp{1}*Z1*T_tmp{1}';
X_b = T_tmp{2}*Z2*T_tmp{2}';

%%%%%%
% Extra
% x_a(ia{1}) = x_b(ia{1});
% x_b(ia{2}) = x_a(ia{2});
% 
% x1(ia{1}) = x2(ia{1});
% x2(ia{2}) = x1(ia{2});
% x1;
% x2;
% 
% X_a(ia{1},ia{1}) = X_b(ia{1},ia{1});
% X_b(ia{2},ia{2}) = X_a(ia{2},ia{2});
%%%%%%

if type == 0
    Zf      = X_a + X_b;
    zf      = x_a + x_b;
    0;
elseif type == 1
    % CI
    ZA      = X_a;
    ZB      = X_b;
    f       = @(w) trace( pinv(w*ZA + (1-w)*ZB) ); % arg (min -f) = arg (max f)
%     omega   = fminbnd(f,0,1,optimset('Display','off'));
    omega   = 0.2;
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
% %     ff      = @(w) trace(pinv(ZA + ZB - ZA*pinv(w*ZB + (1-w)*ZA)*ZB)); % min -f = max f
%     ff      = @(w) trace(-(ZA + ZB - ZA*pinv(w*ZB + (1-w)*ZA)*ZB)); % min -f = max f
%     omega   = fminbnd(ff,0,1,optimset('Display','off'));
    omega   = 0.5;
    Xij     = ZA*inv(omega*ZB + (1-omega)*ZA)*ZB;
    Zf      = X_a + X_b - Xij;

    K = ( X_a - (omega)*Xij )*inv(Zf);
    L = ( X_b - (1 - omega)*Xij )*inv(Zf);
    zf = K*x_a + L*x_b;
    3;
end
xf  = x_tmp;
end