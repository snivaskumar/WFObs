function [xe,Ce] = fuze(z,Z,x,hr,n,Zkk1,xkk1,x_est,typeCZ,typeIFAC,typeWeight,weight);
% function [xe,Ce] = fuze(z,Z,x,hr,n,x_est,typeCZ);
% [xe,Ce] = fuze(z,Z,x,hr,n,x_est,typeCZ);
% typeCZ  1 if Z = Co-Variance,  
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

if strcmp(typeCZ,'C')
    C = Z;
    if typeIFAC == 2
        for i = 1:hr
            Zkk1{i} = inv(Zkk1{i});
        end
    end
elseif strcmp(typeCZ,'Z')
    for i = 1:hr
        C{i} = pinv(Z{i});
    end
end

Zp = zeros(n,n);
for i = 1:hr
    l                   = length(x{i});
    nn                  = nn + l;
    
    if typeIFAC == 2
        Zp(x{i},x{i})       = Zkk1{i};
    end
    if strcmp(typeWeight,'CONSTANT')
%         w = l;
%         w = 1/weight;
        w = weight;
        if i ~= 1
            w = (1-w)/(hr-1);
        end
        Ze(k+1:nn,k+1:nn)   = pinv( (1/w).*C{i} );
    elseif strcmp(typeWeight,'OPTIMAL')
        Ze(k+1:nn,k+1:nn)   = C{i};
%         f = @(w) trace((1/w).*C{i});
%         w = fminbnd(f,0,1,optimset('Display','off'));
%         Ze(k+1:nn,k+1:nn)   = inv( (1/w).*C{i} );
    end
    if strcmp(typeCZ,'C')
        xx(k+1:nn,:)        = z{i};
    elseif strcmp(typeCZ,'Z')
        xx(k+1:nn,:)        = C{i}*z{i};
    end

    k                   = k + l;
end
if strcmp(typeWeight,'OPTIMAL')
    f   = @(w) trace((1/w).*Ze);
    w   = fminbnd(f,0,1,optimset('Display','off'));
    Ze  = inv( (1/w).*Ze );
end

if typeIFAC == 2
    Zp = Zp(x_est,x_est);
    if nargin >= 8
        H = H(:,x_est);
    end
else
    if nargin >= 6
        H = H(:,x_est);
    end
end
H = sparse(H);

% size(Zp)
% size(H'*Ze*H)
if typeIFAC == 2
    Ce          = inv(Zp + H'*Ze*H);
    l_est       = length(x_est);
    [rze,cze]   = size(Ze);
    H           = [eye(l_est,l_est);H];
    Ze          = [Zp, zeros(l_est,cze); zeros(rze,l_est), Ze]; 
    xx          = [xkk1(x_est); xx]; 
else
    Ce = pinv(H'*Ze*H);
%     xe = Ce*H'*Ze*xx;
end
% size(Ze)
% size(H)
% size(xx)
% size(Ce)
xe = Ce*H'*Ze*xx;
end