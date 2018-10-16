% % % clear all
% % % close all
% % % clc
% % % % a = rand(5,1);
% % % % b = a;
% % % % P = rand(5,5);
% % % % P = (P*P')/norm(P*P');
% % % a = [0.4480; 0.5611; 0.0698; 0.4981; 0.7489];
% % % P = [0.3836    0.2449    0.1045    0.3815    0.1293;
% % % 0.2449    0.1616    0.0665    0.2476    0.0809;
% % % 0.1045    0.0665    0.0324    0.1066    0.0317;
% % % 0.3815    0.2476    0.1066    0.3949    0.1252;
% % % 0.1293    0.0809    0.0317    0.1252    0.0556];
% % % Z = inv(P);
% % % z = Z*a;
% % % x = cell(3,1);
% % % x{1} = [1,2,3]';
% % % x{2} = [2,3,4]';
% % % x{3} = [4,5]';
% % % a1 = a(x{1});a2 = a(x{2});a3 = a(x{3});
% % % z1 = z(x{1});z2 = z(x{2});z3 = z(x{3});
% % % P1 = P(x{1},x{1});P2 = P(x{2},x{2});P3 = P(x{3},x{3});
% % % Z1 = Z(x{1},x{1});Z2 = Z(x{2},x{2});Z3 = Z(x{3},x{3});
% % % [zf, Zf] = d_fuze2(a1,a2,P1,P2,x{1},x{2},'CI');
% % % zf{1}
% % % zf{2}
% % % zz{1} = a1;zz{2}=a2;ZZ{1}=P1;ZZ{2}=P2;
% % % [xe,Ce] = d_fuzee(zz,ZZ,x,2,2,0,0,[1:4]','C',1,'OPTIMAL');
% % % xe
% % % Ce

function [xe,Ce] = d_fuze(z,Z,x,hr,sys,Zkk1,xkk1,x_est,typeCZ,typeIFAC,typeWeight);
% [xe,Ce] = d_fuze(z,Z,x,hr,sys,Zkk1,xkk1,x_est,typeCZ,typeIFAC,typeWeight);
% typeCZ  1 if Z = Co-Variance,  
%       2 if Z = Information,

% sys = 1;
n = length(x{sys});
n_sys = length(x{sys});
H = [];
k = 0;
nn = 0;
for i = 1:hr
    if i == sys 
        l               = length(x{i});
        nn              = nn + l;
        H(k+1:nn,:)     = eye(l,l);
        k               = k + l;
    else
        [x12,x1,x2]     = intersect( x{sys},x{i} );
        x12             = sort(x12);
        l               = length( x12 );
        n               = n_sys;
        nn              = nn + l;
        H(k+1:nn,:)     = zeros(l,n);
        H(k+1:nn,x1)  = eye(l,l);
        k               = k + l;
    end
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
        C{i} = inv(Z{i});
    end
end

for i = 1:hr
    [x12,x1,x2]     = intersect( x{sys},x{i} );
    l               = length( x12 );
    nn              = nn + l;
end    
Ze = zeros(nn,nn);
Zp = zeros(n,n);
nn = 0;
for i = 1:hr
    if i == sys
        l                   = length(x{i});
        nn                  = nn + l;

        if typeIFAC == 2
            Zp(x{i},x{i})       = Zkk1{i};
        end
        if strcmp(typeWeight,'CONSTANT')
            w = 1/hr;
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
    else
        [x12,x1,x2]      = intersect( x{sys},x{i} );
        x12              = sort(x12); 
        x1              = sort(x1);
        x2              = sort(x2);
        l               = length( x12 );
        n               = n_sys;
        nn              = nn + l;
        
        if typeIFAC == 2
            Zp(x12,x12)       = Zkk1{i}(x2,x2);
        end
        if strcmp(typeWeight,'CONSTANT')
            w = 1/hr;
            Ze(k+1:nn,k+1:nn)   = inv( (1/w).*C{i}(x2,x2) );
        elseif strcmp(typeWeight,'OPTIMAL')
            Ze(k+1:nn,k+1:nn)   = C{i}(x2,x2);
% % % %     %         f = @(w) trace((1/w).*C{i});
% % % %     %         w = fminbnd(f,0,1,optimset('Display','off'));
% % % %     %         Ze(k+1:nn,k+1:nn)   = inv( (1/w).*C{i} );
        end
        if strcmp(typeCZ,'C')
            xx(k+1:nn,:)        = z{i}(x2);
        elseif strcmp(typeCZ,'Z')
            xx(k+1:nn,:)        = C{i}(x2,x2)*z{i}(x2);
        end
        
        k                   = k + l;
    end
end
if strcmp(typeWeight,'OPTIMAL')
    f   = @(w) trace((w).*Ze);
    w   = fminbnd(f,0,1,optimset('Display','off'));
    Ze  = inv( (1/w).*Ze );
end

if typeIFAC == 2
    Zp = Zp(x_est,x_est);
%     if nargin >= 8
%         H = H(:,x_est);
%     end
else
%     if nargin >= 6
%         H = H(:,x_est);
%     end
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