function [x,d,p, F,D,G,H,Q,R, tur,n, x_est,x_unest, P_unest] = subsystem_turbine(p,Fk,Bk,C,QQ,RR, tur,state,turbLocArray, Subsys_length,RD, Sk1k1);
% [x,d, F,D,G,H,Q,R,y, tur,n, x_est,x_unest, P_unest] = subsystem_turbine(Fk,Bk,C,QQ,RR,yy, tur,stateLocArray,turbLocArray, Sk1k1);

% Model Decomposition
% tic
% aa  = double(abs(Fk)>1e-4);
% p   = symrcm(aa);
% A   = Fk(p,p);
% Bk  = Bk(p,:);
% Ck  = C(:,p);
% state = stateLocArray(p,:);
A       = Fk;
Bk      = Bk;
Ck      = C;
state   = state;
% toc
n = length(Fk);

turbine = turbLocArray;
d = cell(tur,1);
% tic
parfor j = 1:tur
    for i = 1:n
        d{j}(i) = sqrt( (state(i,1)-turbine(j,1))^2 + (state(i,2)-turbine(j,2))^2 );
    end
end
% toc
x = cell(tur,1);
x_ha = cell(tur,1);
% pp = cell(tur,1);

RD;
if Subsys_length <= 5
    Subsys_length = Subsys_length*RD;
else
    Subsys_length = Subsys_length;
end
Subsys_length;

% tic
parfor j = 1:tur
    for i = 1:n
        if d{j}(i)<= (Subsys_length)
            x{j} = [x{j},i];
            x_ha{j} = [x_ha{j};state(i,:)];
%             pp{j} = [pp{j},p(i)];
        end
    end
    x{j} = x{j}';
end
% toc

x_est = x{1};
for i = 2:tur
    x_est = union( x_est, x{i} );
end
x_unest = setdiff([1:n],x_est)';
size(A);
length(x_est);
length(x_unest);

d = cell(tur,1);
% tic
for i = 1:tur
    clear rtmp ctmp tmp
    tmp         = setdiff(x_est,x{i});
    [rtmp ctmp] = find( A(x{i},tmp) );
    ctmp        = unique(ctmp);
    d{i}        = [d{i}, tmp(ctmp)];
end
% toc
for i = 1:tur
    d{i} = unique(d{i});
end

%% Sub-Systems

F   = cell(tur,1);                 % Local A
H   = cell(tur,1);                 % Local C
D   = cell(tur,1);                 % Local Internal Input 
G   = cell(tur,1);                 % Local External Input
Q   = cell(tur,1);
R   = cell(tur,1);
y   = cell(tur,1);
% tic
for i = 1:tur
    F{i}    = A( x{i},x{i} );
    D{i}    = A( x{i},d{i} );
    G{i}    = Bk( x{i},: );
    H{i}    = Ck( :,x{i} );
    R{i}    = RR( :,: );
%     y{i}    = yy( : );
    Q{i}    = QQ( x{i},x{i} );
end
% toc
P_unest = [];
% P_unest = A( x_unest,x_unest )*Sk1k1( x_unest,x_unest )*A( x_unest,x_unest )' ...
%                              + A( x_unest,x_unest )*Sk1k1( x_unest,x_est )*A( x_unest,x_est )' ...
%                              + (A( x_unest,x_unest )*Sk1k1( x_unest,x_est )*A( x_unest,x_est )')' ...
%                              + A( x_unest,x_est )*Sk1k1( x_est,x_est )*A( x_unest,x_est )' ...
%                              + QQ( x_unest,x_unest );

end