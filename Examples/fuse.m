function [zf, X] = fuse(xlkk,Slkk,x,hr);
% [zf, X] = fuse(xlkk,Slkk,x,hr);
% xlkk      Local State vector to be fused
% Slkk      Local Co-Variance Matrix to be fused
% x         Local State vector index
% hr        Number of Local Nodes

% % y_a       = [ 3.312313, 10.24324 ]';
% % y_b       = [ 9.986453, 5.324243, 2.546464 ]';
% 
% % y_b       = [ 9.986453, 5.324243]';
% % ya = [1.13223 10.23131]';
% % yb = [9.96432 5.325672]';
% 
% clear all
% close all
% clc
% 
% hr = 3;
% xlkk        = cell(hr,1);
% xlkk{1}     = [1.13223 10.23131]';
% xlkk{2}     = [9.96432 5.325672]';
% xlkk{3}     = [5.324243 2.546464]';
% Slkk        = cell(hr,1);
% lx          = zeros(hr,1);
% for i = 1:hr
%     lx(i)       = length(xlkk{i});
%     Slkk{i}     = bandmatrix(lx(i),lx(i),lx(i));
% end
% 
% x       = cell(hr,1);
% x{1}    = [1,2]';
% x{2}    = [2,3]';
% x{3}    = [3,4]';

%% 

tic
xf      = x{1};
T_tmp   = cell(hr,1);
cc      = cell(hr,1);
ia      = cell(hr,1);
zf      = xlkk{1};
X       = Slkk{1};
for i = 1:hr-1
x_tmp                           = union( xf,x{i+1} );
[cc{i},ia{i}]                   = setdiff( x_tmp,xf );
T_tmp{i}                        = zeros(max(x_tmp),max(x_tmp));
T_tmp{i}(x_tmp,x_tmp)           = eye(length(x_tmp),length(x_tmp));
T_tmp{i}(cc{i},cc{i})           = zeros(length(cc{i}),length(cc{i}));
% T_tmp{i}                        = T_tmp{i}(x_tmp,setdiff(x_tmp,ia{i}));
T_tmp{i}                        = T_tmp{i}(x_tmp,xf);
T_tmp{i}                        = sparse(T_tmp{i});
x_a                             = T_tmp{i}*zf;

% clear cc ia
[cc{i+1},ia{i+1}]               = setdiff( x_tmp,x{i+1} );
T_tmp{i+1}                      = zeros(max(x_tmp),max(x_tmp));
T_tmp{i+1}(x_tmp,x_tmp)         = eye(length(x_tmp),length(x_tmp));
T_tmp{i+1}(cc{i+1},cc{i+1})     = zeros(length(cc{i+1}),length(cc{i+1}));
% T_tmp{i+1}                      = T_tmp{i+1}(x_tmp,setdiff(x_tmp,ia{i+1}));
T_tmp{i+1}                      = T_tmp{i+1}(x_tmp,x{i+1});
T_tmp{i+1}                        = sparse(T_tmp{i+1});
x_b                             = T_tmp{i+1}*xlkk{i+1};

x_a(ia{i}) = x_b(ia{i});
x_b(ia{i+1}) = x_a(ia{i+1});

[Sa, Da]            = eig(T_tmp{i}*X*T_tmp{i}');
[Sb, Db]            = eig(T_tmp{i+1}*Slkk{i+1}*T_tmp{i+1}');

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

% [Si,Di] = eig(pinv(X_a));
% [Sj,Dj] = eig(pinv(Di^0.5)*pinv(Si)*pinv(X_b)*Si*pinv(Di^0.5));

Dij = zeros(size(Dj));
for ii = 1: length(Dij)
    Dij(ii,ii) = max(1,Dj(ii,ii));
end
% Dij = pinv(Dij);
Xij = Si*(Di^0.5)*Sj*Dij*pinv(Sj)*(Di^0.5)*pinv(Si);

eps = 1e-4;
for ii = 1: length(Dij)
    if abs(1-Dj(ii,ii)) > 10*eps
        cx(ii,ii) = 0; 
    else
        cx(ii,ii) = eps;
    end
end
xij = pinv(pinv(X_a) + pinv(X_b) - 2*pinv(Xij) + 2*cx*eye(size(x_tmp)))*((pinv(X_b) - pinv(Xij) + cx*eye(size(x_tmp)))*x_a + (pinv(X_a) - pinv(Xij) + cx*eye(size(x_tmp)))*x_b);

X   = pinv(pinv(X_a) + pinv(X_b) - pinv(Xij));
zf  = X*(pinv(X_a)*x_a + pinv(X_b)*x_b - pinv(Xij)*xij);

% xij = pinv( X_a + X_b - 2*Xij)*((X_b - Xij)*x_a + (X_a - Xij)*x_b);
% 
% X   = X_a + X_b - Xij;
% zf  = x_a + x_b;

xf  = x_tmp;
clear T_tmp{i} T_tmp{i+1}
end

% xlkk{1}
% xlkk{2}
% xlkk{3}
% zf

toc