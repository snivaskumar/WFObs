function [x,d,p,pp ,F,D,G,H,Q,R,y,hr,n,x_est,x_unest, lc] = subsystem_output( Ak,Bk,Ck,Dk,QQ,RR,yy);
% [x,d,F,D,G,H,Q,R] = subsystem_output( Ak,Bk,Ck,Dk);

% Ak = [11 12 0 0 0; 21 22 0 24 0; 31 0 33 0 0; 0 0 43 0 45; 0 0 0 54 55]
% % Ck = [11 12 13 0 0; 0 22 23 24 0; 0 0 0 34 35]
% Ck = [11 0 0 0 0; 21 0 23 0 0; 0 0 0 0 35]
% Bk = [0 0; 0 0; 0 1; 0 0; 1 0];
% % % % L = bandwidth(A);
% % % % Sk1k1 = bandmatrix(5,L,L);
% % % % QQ = rand(5,5);
% % % % RR = rand(3,3);
% % % % xkk1 = rand(5,1);
% 
% load('/Users/Nivas_Kumar/Desktop/2turb_C_Fk.mat');
% Ck = strucObs.Htt;
% 
% aa = double(abs(Fk)>1e-4);
% spy(aa)
% p = symrcm(aa);
% A = Fk(p,p);
% B = Bk(p,:);
% C = Ck(:,p);
% figure,spy(aa(p,p))
% L = bandwidth(aa(p,p));
% 
% Sk1k1   = bandmatrix(length(p),L,L);
% xkk1 = rand(length(p),1);
% 
% y       = measuredData.sol(strucObs.obs_array);
% % % % Sk1k1   = strucObs.Pk(p,p);      % S(k-1/k-1)
% RR      = strucObs.R_k*randn(strucObs.nrobs,strucObs.nrobs);
% QQ      = strucObs.Q_k(p,p);
% % % % xk1k1   = xk1k1(p);
% % % % xkk1    = solf.x(p);
% % % 
% % % % A = Fk;
% % % % B = Bk;
% % % % C = strucObs.Htt;
% % % % y       = measuredData.sol(strucObs.obs_array);
% % % % Sk1k1   = strucObs.Pk;      % S(k-1/k-1)
% % % % RR      = strucObs.R_k*randn(strucObs.nrobs,strucObs.nrobs);
% % % % QQ      = strucObs.Q_k;
% % % % xkk1    = solf.x;
% 
% [l,u]=bandwidth(double(abs(A) > 1e-4));

aa = double(abs(Ak)>1e-4);
p = symrcm(aa);
A = Ak(p,p);
B = Bk(p,:);
C = Ck(:,p);
D = Dk;

% A = Ak;
% B = Bk;
% C = Ck;
% D = Dk;

n = length(A);
L = bandwidth(A);

[hr hc] = size(C);
[r c] = find(C);
s = unique(c);

%% Localization

if length(s) < n
    xr = [];
    flag = 1;
else
    xr = s;
    flag = 0;
end
k = [];
kk = 1;
nn = 1;
while length(unique(xr')) < n
    for j = 1:length(c)
        xx = find(abs(A(c(j),:))>1e-4);
        xx = [xx, c(j)];
        xx = unique(xx');
        xr = [xr, xx'];
        sxx = length(xx);
        k(kk:kk+sxx-1) = r(j)*ones(1,sxx);
        kk = kk+sxx;
    end
        break;
end

% if strucObs.stateEst || strucObs.measFlow
%     stateLocArray = zeros(strucObs.size_output,2);
%     for iii = 1:strucObs.size_output
%         [~,loci,~]           = WFObs_s_sensors_nr2grid(iii,Wp.mesh);
%         stateLocArray(iii,:) = [loci.x, loci.y];
%     end
% end
% 
% if strucObs.tune.est || strucObs.measPw
%     turbLocArray = zeros(Wp.turbine.N,2);
%     for iii = 1:Wp.turbine.N
%         turbLocArray(iii,:) = [Wp.turbine.Crx(iii),Wp.turbine.Cry(iii)];
%     end
% end
%     
% outputLocArray = [];
% if strucObs.measFlow
% outputLocArray = [outputLocArray; stateLocArray(strucObs.obs_array,:)];
% end
% if strucObs.measPw
% outputLocArray = [outputLocArray; turbLocArray];
% end
% 
% strucObs.r_infl = 1.025;     % Covariance inflation factor (typically 1.00-1.20, no inflation: 1)
% strucObs.f_locl = 'gaspari'; % Localization method: 'off', 'gaspari' (Gaspari-Cohn 1999) or 'heaviside' (Heaviside step function: 0s or 1s)
% strucObs.l_locl = 131;       % Gaspari-Cohn: typically sqrt(10/3)*L with L the cut-off length. Heaviside: cut-off length (m).
% 
% lop     = length(strucObs.obs_array);
% n
% if strucObs.stateEst
%     rho_locl.cross = sparse(lop,n);
%     for iii = 1:strucObs.size_output % Loop over all default states
%         loc1 = stateLocArray(iii,:);
%         for jjj = 1:lop % Loop over all measurements
%             loc2 = outputLocArray(jjj,:);
%             dx = sqrt(sum((loc1-loc2).^2)); % displacement between state and output
%             rho_locl.cross(iii,jjj) = localizationGain( dx, strucObs.f_locl, strucObs.l_locl );
%         end
%     end
%     clear iii jjj dx loc1 loc2
% else
%     rho_locl.cross = [];
% end
% cross = zeros(2137,2137);
% cross(:,1:hr) = rho_locl.cross;
% figure,spy(cross'),xlim([0 hc]),ylim([1 hr])
%% Local States
x = cell(hr,1);
% % % % % if length(s) == n
% % % % %     for i = 1:length(r)
% % % % %         x{r(i)} = [x{r(i)}, c(i)];
% % % % %     end
% % % % % else
% % % % %     for i = 1:length(k)
% % % % %         x{k(i)} = [x{k(i)}, xr(i)];
% % % % %     end
% % % % % end

if flag == 0
    for i = 1:length(r)
        x{r(i)} = [x{r(i)}, c(i)];
    end
else
    for i = 1:length(k)
        x{k(i)} = [x{k(i)}, xr(i)];
    end
end
for i = 1:hr
    x{i} = x{i}';
end

% for i = 1:hr
% %     x{i}    = find(rho_locl.cross(:,i));
%     [rx,cx] = find(p == x{i});
%     [rx,I]  = sort(rx);
%     cx      = cx(I);
%     x{i}    = cx;
%     clear rx cx I
% end

for i = 1:hr
    x{i} = unique(x{i});
end

pp = cell(hr,1);
parfor i = 1:hr
    for j = 1:length(x{i});
        if sum( double( x{j} == p ) ) == 1
            pp{i} = [pp{i};x{j}];
        end
    end
end

x_est = x{1};
for i = 2:hr
    x_est = union( x_est, x{i} );
end
x_unest = setdiff([1:n],x_est)';
size(A);
length(x_est);
length(x_unest);
%% Internal Input
d = cell(hr,1);

% for i = 1:hr
%     dA = [];
%     for j = 1:n
%         if sum(j == x{i}) == 0
%             dA = [dA A(x{i},j)];
%             if sum(A(x{i},j)~=0) > 0
%                 d{i} = [d{i}, j];
%             end
%         end
%     end
%     dA;
% end

% Same as the above one (but with inbuilt codes)
ddd = zeros(n,1);
for i = 1:hr
    clear rtmp ctmp tmp
    tmp         = setdiff(x_est,x{i});
%     tmp         = find(~rho_locl.cross(:,i));
    [rtmp ctmp] = find( A(x{i},tmp) );
    ctmp        = unique(ctmp);
    d{i}        = [d{i}, tmp(ctmp)];
%     ddd         = [ddd, d{i}];
end

for i = 1:hr
    d{i} = unique(d{i});
end

%% Sub-Systems

F   = cell(hr,1);                 % Local A
Fl  = cell(hr,1);
H   = cell(hr,1);                 % Local C
D   = cell(hr,1);                 % Local Internal Input 
G   = cell(hr,1);                 % Local External Input
T   = cell(hr,1);
Q   = cell(hr,1);
R   = cell(hr,1);
y   = cell(hr,1);
for i = 1:hr
%     Ttmp	= eye(n,n);
%     Ttmp(setdiff([1:n],x{i}),setdiff([1:n],x{i})) = 0;
%     T{i}    = Ttmp;
%     T{i}    = sparse(T{i});
    F{i}    = A( x{i},x{i} );
%     Fl{i}   = T{i}*A;
%     Fl{i}   = sparse(Fl{i});
    D{i}    = A( x{i},d{i} );
    G{i}    = B( x{i},: );
    H{i}    = C( i,x{i} );
    y{i}    = yy(i);
    Q{i}    = QQ( x{i},x{i} );
    R{i}    = RR( i,i );
end

end