clear all
close all
clc

load('/Users/Nivas_Kumar/Desktop/2turb_C_Fk.mat');
C = strucObs.Htt;
aa = double(abs(Fk)>1e-4);
p = symrcm(aa);
A = Fk(p,p);
Bk = Bk(p,:);
Ck = C(:,p);

[i,j] = find(Fk);
bw = max(i-j) + 1
[i,j] = find(A);
bw = max(i-j) + 1

state = stateLocArray(p,:);
% state = stateLocArray;
turbine = turbLocArray;
turb1 = turbine(1,:);
turb2 = turbine(2,:);
d = sqrt( (turb2(1)-turb1(1))^2 + (turb2(2)-turb1(2))^2 )

for i = 1:2137
    d1(i) = sqrt( (state(i,1)-turb1(1))^2 + (state(i,2)-turb1(2))^2 );
    d2(i) = sqrt( (state(i,1)-turb2(1))^2 + (state(i,2)-turb2(2))^2 );
end
% plot(d1),hold on, plot(d2), hold off
x_ho = cell(2,1);
x_ha = cell(2,1);
pp = cell(2,1);
pptmp = [];
for i = 1:2137
    if d1(i)<= (d/1)
        x_ho{1} = [x_ho{1},i];
        x_ha{1} = [x_ha{1};state(i,:)];
        pp{1} = [pp{1},p(i)];
    end
    if d2(i)<= (d/1)
        x_ho{2} = [x_ho{2},i];
        x_ha{2} = [x_ha{2};state(i,:)];
        pp{2} = [pp{2},p(i)];
    end
    if ( d1(i) > (d/1) ) && ( d2(i) > (d/1) )
        pptmp = [pptmp;state(i,:)];
    end
end
% figure, plot(x_ho{1},x_ho{1},'o'), hold on, plot(x_ho{2},x_ho{2},'+'), hold off

figure, plot(state(:,2),state(:,1),'+')
hold on, plot(turb1(2),turb1(1),'s','LineWidth',2)
hold on, plot(turb2(2),turb2(1),'s','LineWidth',2), hold off
xlabel('y-direction'),ylabel('x-direction');

figure, plot(x_ha{1}(:,2),x_ha{1}(:,1),'+')
hold on, plot(turb1(2),turb1(1),'s','LineWidth',2)
hold on, plot(turb2(2),turb2(1),'s','LineWidth',2)
hold on, plot(x_ha{2}(:,2),x_ha{2}(:,1),'s')
hold on, plot(pptmp(:,2),pptmp(:,1),'*'), hold off
xlabel('y-direction'),ylabel('x-direction');
legend('States of subsystem-1','Turbine 1','Turbine 2','States of subsystem-2','Unestimated States')
title('Flow Field');

%%

C = strucObs.Htt;
tur         = 2;

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
x = sort(x)

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

lc = cell(tur,1);
for i = 1:tur
    lc{i} = [1:l]'.*qq{i};
    lc{i} = find(lc{i});
    lc{i};
    Ck(lc{i},:);
end

