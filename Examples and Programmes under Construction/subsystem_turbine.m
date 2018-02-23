function [x,d, F,D,G,H,Q,R, tur,n, x_est,x_unest] = subsystem_turbine(Fk,Bk,C,QQ,RR, tur,stateLocArray);
% [x,d,F,D,G,H,Q,R,hr,n,x_est,x_unest] = subsystem_turbine(Fk,Bk,C,QQ,RR, tur,stateLocArray);

% load('/Users/Nivas_Kumar/Desktop/2turb_C_Fk.mat');
% C = strucObs.Htt;
% tur         = 2;

[~,loc_C,~]     = find(C);
lock            = stateLocArray(loc_C,1);

aa = double(abs(Fk)>1e-4);
p = symrcm(aa);
A = Fk(p,p);
Bk = Bk(p,:);
Ck = C(:,p);

[~,loc_Ck,~]    = find(Ck);
[Y,I] = sort( p(loc_Ck) );
stateLocArray   = stateLocArray(p,:);
loc             = stateLocArray(loc_Ck,1);
loc = loc(I);

x           = [];
k           = 0;
kk          = 1; 
xx          = loc;
kk          = length(x);
while kk ~=tur
    l           = length(xx);
        
    max_loc     = max(xx);
    min_loc     = min(xx);
        
    max_q = numel(num2str(round(max_loc)));
    min_q = numel(num2str(round(min_loc)));
        
    max_q = 10^(max_q - 1);
    min_q = 10^(min_q - 1);
        
    max_con_max = max_loc + max_loc/max_q;
    max_con_min = max_loc - max_loc/max_q;
        
    min_con_max = min_loc + min_loc/min_q;
    min_con_min = min_loc - min_loc/min_q;
        
    x(kk + 1)   = max_loc;
    if (max_con_min < min_loc) && (max_loc < min_con_max);
    else
        x(kk + 2)   = min_loc;
    end
    xx1 = [];
    for ii = 1:l
        if (min_con_min < xx(ii)) && (xx(ii) < min_con_max)
        else
        	if (max_con_min < xx(ii)) && (xx(ii) < max_con_max)
            else
            	xx1 = [xx1,xx(ii)];
            end
        end
    end
    xx = xx1;
    k  = k + 1;
    kk = length(x);
end
x = sort(x);

loc_tmp = [];
loc_x = cell(tur,1);
for i = 1:tur
    if i == 1
        avg         = (x(i + 1) - x(i))/2;
        x_con_max   = x(i) + avg;
        q           = (loc <= x_con_max);      
    elseif (i~=1) && (i~=tur)
        avg1        = (x(i) - x(i-1))/2;
        x_con_min   = x(i) - avg1;
        avg2        = (x(i + 1) - x(i))/2;
        x_con_max   = x(i) + avg2;
        q           = (loc <= x_con_max)&&(x_con_min < loc);    
    else
        avg         = (x(i) - x(i-1))/2;
        x_con_min   = x(i) - avg;
        q           = (loc > x_con_min);    
    end
    qq{i} = double(q);
    [r c]   = find(Ck(q,:));
    s{i}    = unique(c);
    loc_tmp = union(loc_tmp,s{i});
end

n = length(A);
l = length(loc_tmp);
x = cell(tur,1);
if l < n
    xrr = cell(tur,1);
    for i = 1:tur
        xxx = [];
        for j = 1:length(s{i})
    	[q,w] = find(abs(A(s{i}(j),:))>1e-4);
        xx = [w'; s{i}(j)];
        x{i} = [x{i}; xx];
        x{i} = unique(x{i});
        end
    end
else
    x = s;
end

d = cell(tur,1);
for i = 1:tur
    clear rtmp ctmp tmp
    tmp         = setdiff([1:n],x{i});
    [rtmp ctmp] = find( A(x{i},tmp) );
    ctmp        = unique(ctmp);
    d{i}        = [d{i}, tmp(ctmp)];
end
for i = 1:tur
    d{i} = unique(d{i});
end

x_est = x{1};
for i = 2:tur
    x_est = union( x_est, x{i} );
end
x_unest = setdiff([1:n],x_est)';
size(A);
length(x_est);
length(x_unest);

%% Sub-Systems

F   = cell(tur,1);                 % Local A
Fl  = cell(tur,1);
H   = cell(tur,1);                 % Local C
D   = cell(tur,1);                 % Local Internal Input 
G   = cell(tur,1);                 % Local External Input
T   = cell(tur,1);
Q   = cell(tur,1);
R   = cell(tur,1);
for i = 1:tur
%     Ttmp	= eye(n,n);
%     Ttmp(setdiff([1:n],x{i}),setdiff([1:n],x{i})) = 0;
%     T{i}    = Ttmp;
%     T{i}    = sparse(T{i});
    F{i}    = A( x{i},x{i} );
%     Fl{i}   = T{i}*A;
%     Fl{i}   = sparse(Fl{i});
    D{i}    = A( x{i},d{i} );
    G{i}    = Bk( x{i},: );
% % %     H{i}    = C( i,x{i} );
    H{i}    = Ck( i,: );
    Q{i}    = QQ( x{i},x{i} );
    R{i}    = RR( i,i );
end

end