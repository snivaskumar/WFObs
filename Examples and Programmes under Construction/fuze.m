function [xe,Ce] = fuze(z,Z,x,hr,n,x_est);
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
for i = 1:hr
    C{i} = inv(Z{i});
end

for i = 1:hr
    l                   = length(x{i});
    nn                  = nn + l;
    
    Ze(k+1:nn,k+1:nn)   = inv( (1/l).*C{i} );
    xx(k+1:nn,:)        = C{i}*z{i};

    k                   = k + l;
end

if nargin == 6
    H = H(:,x_est);
end

Ce = inv(H'*Ze*H);
xe = Ce*H'*Ze*xx;

end