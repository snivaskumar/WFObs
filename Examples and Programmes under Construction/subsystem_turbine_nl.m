function [x,d, tur,n, x_est,x_unest] = subsystem_turbine_nl(tur,stateLocArray,turbLocArray,Subsys_length);

% Model Decomposition
state = stateLocArray;
n = length( state );

turbine = turbLocArray;
turb1 = turbine(1,:);
turb2 = turbine(2,:);
d = sqrt( (turb2(1)-turb1(1))^2 + (turb2(2)-turb1(2))^2 );

for i = 1:n
    d1(i) = sqrt( (state(i,1)-turb1(1))^2 + (state(i,2)-turb1(2))^2 );
    d2(i) = sqrt( (state(i,1)-turb2(1))^2 + (state(i,2)-turb2(2))^2 );
end
x = cell(tur,1);
x_ha = cell(tur,1);
% pp = cell(tur,1);

if Subsys_length == 1
    Subsys_length = d;
elseif Subsys_length == 2
    Subsys_length = d/2;
else
    Subsys_length = Subsys_length;
end
for i = 1:n
    if d1(i)<= (Subsys_length)
        x{1} = [x{1},i];
        x_ha{1} = [x_ha{1};state(i,:)];
%         pp{1} = [pp{1},p(i)];
    end
    if d2(i)<= (Subsys_length)
        x{2} = [x{2},i];
        x_ha{2} = [x_ha{2};state(i,:)];
%         pp{2} = [pp{2},p(i)];
    end
end

x{1} = x{1}';
x{2} = x{2}';

x_est = x{1};
for i = 2:tur
    x_est = union( x_est, x{i} );
end
x_unest = setdiff([1:n],x_est)';
% size(A);
length(x_est);
length(x_unest);

d = cell(tur,1);
for i = 1:tur
    d{i} = setdiff(x_est, x{i});
    d{i} = unique(d{i});
end
end