function [x,d,p,pp, F,D,G,H,Q,R,y, tur,n, x_est,x_unest, P_unest] = subsystem_turbine(Fk,Bk,C,QQ,RR,yy, tur,stateLocArray,turbLocArray, Sk1k1);
% [x,d, F,D,G,H,Q,R,y, tur,n, x_est,x_unest, P_unest] = subsystem_turbine(Fk,Bk,C,QQ,RR,yy, tur,stateLocArray,turbLocArray, Sk1k1);

% Model Decomposition

% load('/Users/Nivas_Kumar/Desktop/2turb_C_Fk.mat');
% C = strucObs.Htt;
% tur         = 2;

aa  = double(abs(Fk)>1e-4);
p   = symrcm(aa);
A   = Fk(p,p);
Bk  = Bk(p,:);
Ck  = C(:,p);
state = stateLocArray(p,:);
n = length(Fk);

turbine = turbLocArray;
turb1 = turbine(1,:);
turb2 = turbine(2,:);
d = sqrt( (turb2(1)-turb1(1))^2 + (turb2(2)-turb1(2))^2 );

for i = 1:n
    d1(i) = sqrt( (state(i,1)-turb1(1))^2 + (state(i,2)-turb1(2))^2 );
    d2(i) = sqrt( (state(i,1)-turb2(1))^2 + (state(i,2)-turb2(2))^2 );
end
% plot(d1),hold on, plot(d2), hold off
x = cell(tur,1);
x_ha = cell(tur,1);
pp = cell(tur,1);
for i = 1:n
    if d1(i)<= (d/1)
        x{1} = [x{1},i];
        x_ha{1} = [x_ha{1};state(i,:)];
        pp{1} = [pp{1},p(i)];
    end
    if d2(i)<= (d/1)
        x{2} = [x{2},i];
        x_ha{2} = [x_ha{2};state(i,:)];
        pp{2} = [pp{2},p(i)];
    end
end

x{1} = x{1}';
x{2} = x{2}';

% [~,loc_C,~]     = find(C);
% lock            = stateLocArray(loc_C,1);
% 
% aa = double(abs(Fk)>1e-4);
% p = symrcm(aa);
% A = Fk(p,p);
% Bk = Bk(p,:);
% Ck = C(:,p);
% 
% [~,loc_Ck,~]    = find(Ck);
% [Y,I] = sort( p(loc_Ck) );
% stateLocArray   = stateLocArray(p,:);
% loc             = stateLocArray(loc_Ck,1);
% loc = loc(I);
% 
% x           = [];
% k           = 0;
% kk          = 1; 
% xx          = loc;
% kk          = length(x);
% while kk ~=tur
%     l           = length(xx);
%         
%     max_loc     = max(xx);
%     min_loc     = min(xx);
%         
%     max_q = numel(num2str(round(max_loc)));
%     min_q = numel(num2str(round(min_loc)));
%         
%     max_q = 10^(max_q - 1);
%     min_q = 10^(min_q - 1);
%         
%     max_con_max = max_loc + max_loc/max_q;
%     max_con_min = max_loc - max_loc/max_q;
%         
%     min_con_max = min_loc + min_loc/min_q;
%     min_con_min = min_loc - min_loc/min_q;
%         
%     x(kk + 1)   = max_loc;
%     if (max_con_min < min_loc) && (max_loc < min_con_max);
%     else
%         x(kk + 2)   = min_loc;
%     end
%     xx1 = [];
%     for ii = 1:l
%         if (min_con_min < xx(ii)) && (xx(ii) < min_con_max)
%         else
%         	if (max_con_min < xx(ii)) && (xx(ii) < max_con_max)
%             else
%             	xx1 = [xx1,xx(ii)];
%             end
%         end
%     end
%     xx = xx1;
%     k  = k + 1;
%     kk = length(x);
% end
% x = sort(x);
% 
% loc_tmp = [];
% loc_x = cell(tur,1);
% for i = 1:tur
%     if i == 1
%         avg         = (x(i + 1) - x(i))/2;
%         x_con_max   = x(i) + avg;
%         q           = (loc <= x_con_max);      
%     elseif (i~=1) && (i~=tur)
%         avg1        = (x(i) - x(i-1))/2;
%         x_con_min   = x(i) - avg1;
%         avg2        = (x(i + 1) - x(i))/2;
%         x_con_max   = x(i) + avg2;
%         q           = (loc <= x_con_max)&&(x_con_min < loc);    
%     else
%         avg         = (x(i) - x(i-1))/2;
%         x_con_min   = x(i) - avg;
%         q           = (loc > x_con_min);    
%     end
%     qq{i} = double(q);
%     [r c]   = find(Ck(q,:));
%     s{i}    = unique(c);
%     loc_tmp = union(loc_tmp,s{i});
% end
% 
% n = length(A);
% l = length(loc_tmp);
% x = cell(tur,1);
% if l < n
%     xrr = cell(tur,1);
%     for i = 1:tur
%         xxx = [];
%         for j = 1:length(s{i})
%     	[q,w] = find(abs(A(s{i}(j),:))>1e-4);
%         xx = [w'; s{i}(j)];
%         x{i} = [x{i}; xx];
%         x{i} = unique(x{i});
%         end
%     end
% else
%     x = s;
% end
% 
% lc = cell(tur,1);
% for i = 1:tur
%     lc{i} = [1:l]'.*qq{i};
%     lc{i} = find(lc{i});
% end

x_est = x{1};
for i = 2:tur
    x_est = union( x_est, x{i} );
end
x_unest = setdiff([1:n],x_est)';
size(A);
length(x_est);
length(x_unest);

d = cell(tur,1);
% for i = 1:tur
%     clear rtmp ctmp tmp
%     tmp         = setdiff([1:n],x{i});
%     [rtmp ctmp] = find( A(x{i},tmp) );
%     ctmp        = unique(ctmp);
%     d{i}        = [d{i}, tmp(ctmp)];
% end

for i = 1:tur
    clear rtmp ctmp tmp
    tmp         = setdiff(x_est,x{i});
    [rtmp ctmp] = find( A(x{i},tmp) );
    ctmp        = unique(ctmp);
    d{i}        = [d{i}, tmp(ctmp)];
end
for i = 1:tur
    d{i} = unique(d{i});
end


% % % % for i = 1:tur
% % % %     x{i} = union( x{i}, x_unest );
% % % % end
% % % % x_est = [1:n]';
% % % % x_unest = [];

%% Sub-Systems

F   = cell(tur,1);                 % Local A
H   = cell(tur,1);                 % Local C
D   = cell(tur,1);                 % Local Internal Input 
G   = cell(tur,1);                 % Local External Input
Q   = cell(tur,1);
R   = cell(tur,1);
y   = cell(tur,1);
for i = 1:tur
    F{i}    = A( x{i},x{i} );
    D{i}    = A( x{i},d{i} );
    G{i}    = Bk( x{i},: );
    
%     H{i}    = Ck( i,x{i} );
%     R{i}    = RR( i,i );
%     y{i}    = yy( i );
%     H{i}    = Ck( lc{i},x{i} );
%     R{i}    = RR( lc{i},lc{i} );
%     y{i}    = yy( lc{i} );
    H{i}    = Ck( :,x{i} );
    R{i}    = RR( :,: );
    y{i}    = yy( : );
    
    Q{i}    = QQ( x{i},x{i} );
end
% tmp1 = [1:5,11]';
% tmp2 = [6:10,15]';

%     H{1}    = Ck( tmp1,x{1} );
%     R{1}    = RR( tmp1,tmp1 );
%     y{1}    = yy( tmp1 );
%     
%     H{2}    = Ck( tmp2,x{2} );
%     R{2}    = RR( tmp2,tmp2 );
%     y{2}    = yy( tmp2 );

P_unest = [];
% P_unest = A( x_unest,x_unest )*Sk1k1( x_unest,x_unest )*A( x_unest,x_unest )' ...
%                              + A( x_unest,x_unest )*Sk1k1( x_unest,x_est )*A( x_unest,x_est )' ...
%                              + (A( x_unest,x_unest )*Sk1k1( x_unest,x_est )*A( x_unest,x_est )')' ...
%                              + A( x_unest,x_est )*Sk1k1( x_est,x_est )*A( x_unest,x_est )' ...
%                              + QQ( x_unest,x_unest );

end