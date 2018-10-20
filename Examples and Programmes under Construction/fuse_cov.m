function [X, X_a, X_b, Xij, xf, xab_size, T1, T2, ia1, ia2, cx] = fuse_cov(P1,P2,x1,x2);
% [X, X_a, X_b, Xij, xf, xab_size, T_tmp{1}, T_tmp{2}, ia{1}, ia{2}, cx] = fuse_cov(P1,P2,x1,x2);
% xlkk      Local State vector to be fused
% Slkk      Local Co-Variance Matrix to be fused
% x         Local State vector index
% hr        Number of Local Nodes

xf      = x1;
T_tmp   = cell(2,1);
cc      = cell(2,1);
ia      = cell(2,1);
X       = P1;

x_tmp                           = union( xf,x2 );
xab_size                        = x_tmp;

[cc{1},ia{1}]                   = setdiff( x_tmp,xf );
T_tmp{1}                        = zeros(max(x_tmp),max(x_tmp));
T_tmp{1}(x_tmp,x_tmp)           = eye(length(x_tmp),length(x_tmp));
T_tmp{1}(cc{1},cc{1})           = zeros(length(cc{1}),length(cc{1}));
% T_tmp{i}                        = T_tmp{i}(x_tmp,setdiff(x_tmp,ia{i}));
T_tmp{1}                        = T_tmp{1}(x_tmp,xf);
T_tmp{1}                        = sparse(T_tmp{1});
T1                              = T_tmp{1};
% x_a                             = T_tmp{1}*zf;

% clear cc ia
[cc{2},ia{2}]                   = setdiff( x_tmp,x2 );
T_tmp{2}                        = zeros(max(x_tmp),max(x_tmp));
T_tmp{2}(x_tmp,x_tmp)           = eye(length(x_tmp),length(x_tmp));
T_tmp{2}(cc{2},cc{2})           = zeros(length(cc{2}),length(cc{2}));  
% T_tmp{i+1}                      = T_tmp{i+1}(x_tmp,setdiff(x_tmp,ia{i+1}));
T_tmp{2}                        = T_tmp{2}(x_tmp,x2);
T_tmp{2}                        = sparse(T_tmp{2});
T2                              = T_tmp{2};
% x_b                             = T_tmp{2}*z2;

ia1 = ia{1};
ia2 = ia{2};
% x_a(ia{1}) = x_b(ia{1});
% x_b(ia{2}) = x_a(ia{2});

[Sa, Da]            = eig(T_tmp{1}*X*T_tmp{1}');
[Sb, Db]            = eig(T_tmp{2}*P2*T_tmp{2}');

    for ii = 1:length(x_tmp)
        for jj =1:length(x_tmp)
            if Db(ii,jj) == 0
                Gb(ii,jj) = 1;
            else
                Gb(ii,jj) = Db(ii,jj);
            end
        end
    end

    for ii = 1:length(x_tmp)
        for jj =1:length(x_tmp)
            if Da(ii,jj) == 0
                Ga(ii,jj) = 1;
            else
                Ga(ii,jj) = Da(ii,jj);
            end
        end
    end

X_b = Sa*(Ga^0.5)*Sb*Gb*pinv(Sb)*(Ga^0.5)*pinv(Sa);
X_a = Sb*(Gb^0.5)*Sa*Ga*pinv(Sa)*(Gb^0.5)*pinv(Sb);

[Si,Di] = eig(X_a);
[Sj,Dj] = eig(pinv(Di^0.5)*pinv(Si)*X_b*Si*pinv(Di^0.5));

Dij = zeros(size(Dj));
for ii = 1: length(Dij)
    Dij(ii,ii) = max(1,Dj(ii,ii));
end

Xij = Si*(Di^0.5)*Sj*Dij*pinv(Sj)*(Di^0.5)*pinv(Si);

eps = 1e-4;
for ii = 1: length(Dij)
    if abs(1-Dj(ii,ii)) > 10*eps
        cx(ii,ii) = 0; 
    else
        cx(ii,ii) = eps;
    end
end
% xij = pinv(pinv(X_a) + pinv(X_b) - 2*pinv(Xij) + 2*cx*eye(size(x_tmp)))*((pinv(X_b) - pinv(Xij) + cx*eye(size(x_tmp)))*x_a + (pinv(X_a) - pinv(Xij) + cx*eye(size(x_tmp)))*x_b);

X   = pinv(pinv(X_a) + pinv(X_b) - pinv(Xij));
% zf  = X*(pinv(X_a)*x_a + pinv(X_b)*x_b - pinv(Xij)*xij);
xf  = x_tmp;