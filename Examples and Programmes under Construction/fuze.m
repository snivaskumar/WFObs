function [xe,Ce] = fuze(z,Z,x,hr,n,x_est,type);
% [xe,Ce] = fuze(z,Z,x,hr,n,x_est,type);
% type  1 if Z = Co-Variance,  
%       2 if Z = Information,

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

Ze  = [];
nn  = 0;
k   = 0;  

if type == 1
    C = Z;
else
    for i = 1:hr
        C{i} = inv(Z{i});
    end
end
    

for i = 1:hr
    l                   = length(x{i});
    nn                  = nn + l;

    Ze(k+1:nn,k+1:nn)   = inv( (1/l).*C{i} );
    if type == 1
        xx(k+1:nn,:)        = z{i};
    else
        xx(k+1:nn,:)        = C{i}*z{i};
    end

    k                   = k + l;
end

if nargin >= 6
    H = H(:,x_est);
end

Ce = inv(H'*Ze*H);
xe = Ce*H'*Ze*xx;

end