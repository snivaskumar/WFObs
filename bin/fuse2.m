function [zf, Zf, xf] = fuse2(z1,z2,Z1,Z2,x1,x2,type);
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
    Zf      = pinv(pinv(X_a) + pinv(X_b));
    zf      = Zf*(pinv(X_a)*x_a + pinv(X_b)*x_b);
    0;
elseif type == 1
    % CI
    ZA      = X_a;
    ZB      = X_b;
    f       = @(w) trace( pinv(w*pinv(ZA) + (1-w)*pinv(ZA)) ); % arg (min -f) = arg (max f)
    omega   = fminbnd(f,0,1,optimset('Display','off'));
    Zf      = pinv(omega*pinv(ZA) + (1-omega)*pinv(ZA));
    zf      = Zf*(omega*pinv(ZA)*x_a + (1-omega)*pinv(ZB)*x_b);
    1;
elseif type == 2
    % EI
    [Si,Di] = eig(X_a);
    [Sj,Dj] = eig(pinv(Di^0.5)*pinv(Si)*X_b*Si*pinv(Di^0.5));

    Dij = zeros(size(Dj));
    for ii = 1: length(Dij)
        Dij(ii,ii) = max(1,Dj(ii,ii));
    end
    Xij = Si*(Di^0.5)*Sj*Dij*pinv(Sj)*(Di^0.5)*pinv(Si);

%         w1 = (X_b - 0.5*Xij)*pinv(X_a + X_b - Xij);
%         w2 = (X_a - 0.5*Xij)*pinv(X_a + X_b - Xij);
        
    w1 = pinv(pinv(X_a) + pinv(X_b) - 2*pinv(Xij))*(pinv(X_b) - pinv(Xij));
    w2 = pinv(pinv(X_a) + pinv(X_b) - 2*pinv(Xij))*(pinv(X_a) - pinv(Xij));
    xij = w1*x_a + w2*x_b;

    Zf  = pinv(pinv(X_a) + pinv(X_b) - pinv(Xij));
    zf  = Zf*(inv(X_a)*x_a + inv(X_b)*x_b - inv(X_ij)*xij);
    2;
elseif type == 3    
    % ICI
    ZA      = X_a;
    ZB      = X_b;
    ff      = @(w) trace( (pinv(ZA) + pinv(ZB) - pinv(w*ZA + (1-w)*ZB))); % min -f = max f
    omega   = fminbnd(ff,0,1,optimset('Display','off'));
    Xij     = omega*ZA + (1-omega)*ZB;
    Zf      = pinv(pinv(X_a) + pinv(X_b) - pinv(Xij));

    K = Zf*(pinv(X_a) - omega*pinv(Xij));
    L = Zf*(pinv(X_b) - (1 - omega)*pinv(Xij));
    zf = K*x_a + L*x_b;
    3;
end
xf  = x_tmp;
end


% % function [zf, X, xf] = fuse2(z1,z2,P1,P2,x1,x2);
% % % [zf, X, xf] = fuse2(z1,z2,P1,P2,x1,x2);
% % % xlkk      Local State vector to be fused
% % % Slkk      Local Co-Variance Matrix to be fused
% % % x         Local State vector index
% % % hr        Number of Local Nodes
% % 
% % % tic
% % xf      = x1;
% % T_tmp   = cell(2,1);
% % cc      = cell(2,1);
% % ia      = cell(2,1);
% % zf      = z1;
% % X       = P1;
% % 
% % x_tmp                           = union( xf,x2 );
% % [cc{1},ia{1}]                   = setdiff( x_tmp,xf );
% % T_tmp{1}                        = zeros(max(x_tmp),max(x_tmp));
% % T_tmp{1}(x_tmp,x_tmp)           = eye(length(x_tmp),length(x_tmp));
% % T_tmp{1}(cc{1},cc{1})           = zeros(length(cc{1}),length(cc{1}));
% % % T_tmp{i}                        = T_tmp{i}(x_tmp,setdiff(x_tmp,ia{i}));
% % T_tmp{1}                        = T_tmp{1}(x_tmp,xf);
% % T_tmp{1}                        = sparse(T_tmp{1});
% % x_a                             = T_tmp{1}*zf;
% % 
% % % clear cc ia
% % [cc{2},ia{2}]                   = setdiff( x_tmp,x2 );
% % T_tmp{2}                        = zeros(max(x_tmp),max(x_tmp));
% % T_tmp{2}(x_tmp,x_tmp)           = eye(length(x_tmp),length(x_tmp));
% % T_tmp{2}(cc{2},cc{2})           = zeros(length(cc{2}),length(cc{2}));  
% % % T_tmp{i+1}                      = T_tmp{i+1}(x_tmp,setdiff(x_tmp,ia{i+1}));
% % T_tmp{2}                        = T_tmp{2}(x_tmp,x2);
% % T_tmp{2}                        = sparse(T_tmp{2});
% % x_b                             = T_tmp{2}*z2;
% % 
% % x_a(ia{1}) = x_b(ia{1});
% % x_b(ia{2}) = x_a(ia{2});
% % 
% % [Sa, Da]            = eig(T_tmp{1}*X*T_tmp{1}');
% % [Sb, Db]            = eig(T_tmp{2}*P2*T_tmp{2}');
% % 
% % % Sa = real(Sa);
% % % Da = real(Da);
% % % Sb = real(Sb);
% % % Db = real(Db);
% % 
% %     for ii = 1:length(x_tmp)
% %         for jj =1:length(x_tmp)
% %             if Db(ii,jj) == 0
% %                 Gb(ii,jj) = 1;
% %             else
% %                 Gb(ii,jj) = Db(ii,jj);
% %             end
% %         end
% %     end
% % 
% %     for ii = 1:length(x_tmp)
% %         for jj =1:length(x_tmp)
% %             if Da(ii,jj) == 0
% %                 Ga(ii,jj) = 1;
% %             else
% %                 Ga(ii,jj) = Da(ii,jj);
% %             end
% %         end
% %     end
% % 
% % X_b = Sa*(Ga^0.5)*Sb*Gb*pinv(Sb)*(Ga^0.5)*pinv(Sa);
% % X_a = Sb*(Gb^0.5)*Sa*Ga*pinv(Sa)*(Gb^0.5)*pinv(Sb);
% % 
% % [Si,Di] = eig(X_a);
% % [Sj,Dj] = eig(pinv(Di^0.5)*pinv(Si)*X_b*Si*pinv(Di^0.5));
% % 
% % 
% % % Si = real(Si);
% % % Di = real(Di);
% % % Sj = real(Sj);
% % % Dj = real(Dj);
% % 
% % Dij = zeros(size(Dj));
% % for ii = 1: length(Dij)
% %     Dij(ii,ii) = max(1,Dj(ii,ii));
% % end
% % 
% % Xij = Si*(Di^0.5)*Sj*Dij*pinv(Sj)*(Di^0.5)*pinv(Si);
% % 
% % eps = 1e-4;
% % for ii = 1: length(Dij)
% %     if abs(1-Dj(ii,ii)) > 10*eps
% %         cx(ii,ii) = 0; 
% %     else
% %         cx(ii,ii) = eps;
% %     end
% % end
% % nn = length(x_tmp);
% % xij = pinv(pinv(X_a) + pinv(X_b) - 2*pinv(Xij) + 2*cx*eye(nn,nn))*((pinv(X_b) - pinv(Xij) + cx*eye(nn,nn))*x_a + (pinv(X_a) - pinv(Xij) + cx*eye(nn,nn))*x_b);
% % 
% % X   = pinv(pinv(X_a) + pinv(X_b) - pinv(Xij));
% % zf  = X*(pinv(X_a)*x_a + pinv(X_b)*x_b - pinv(Xij)*xij);
% % xf  = x_tmp;