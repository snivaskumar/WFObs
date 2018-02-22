function [x,C] = fuze(z,Z,x,hr,n,x_est);
H = [];
k = 0;
nn = 0;
for i = 1:hr
    l               = length(x{i});
    nn              = nn + l;
    H(k+1:nn,:)     = zeros(l,n);
    H(k+1:nn,x{i})  = eye(l,l);
    k               = k + l;
end

C = [];
nn = 0;
k = 0;
for i = 1:hr
    l               = length(x{i});
    nn              = nn + l;
%     C(k+1:nn,k+1:nn)= (l).*Z{i};
%     C(k+1:nn,k+1:nn)= (1/l).*Z{i};
%     C(k+1:nn,k+1:nn)= Z{i};
%     xx(k+1:nn,:)    = z{i};

    CC              = inv(Z{i});     
%     C(k+1:nn,k+1:nn)= (1/l).*CC;
    C(k+1:nn,k+1:nn)= inv( (1/l).*CC );
    xx(k+1:nn,:)    = CC*z{i};

    k               = k + l;
end

% Ze = inv(H'*inv(C)*H);
% ze = Ze*H'*inv(C)*xx;

Ze = pinv(H'*(C)*H);
ze = Ze*H'*(C)*xx;

% Ce = inv(Ze);
% xe = Ce*ze;

% C = Ze;
% x = ze;

C = Ze(x_est,x_est);
x = ze(x_est);
end