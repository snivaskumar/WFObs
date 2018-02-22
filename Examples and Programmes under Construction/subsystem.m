function [xkk Pkk] = subsystem( Ak,Bk,Ck,Dk, y, xkk1,xk1k1,Sk1k1, QQ,RR, type );

% function [xkk Pkk] = subsystem( Ak,Bk,Ck,Dk, y, xkk1,xk1k1,Sk1k1, QQ,RR, strucObs,Wp );
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
    tmp         = setdiff([1:n],x{i});
%     tmp         = find(~rho_locl.cross(:,i));
    [rtmp ctmp] = find( A(x{i},tmp) );
    ctmp        = unique(ctmp);
    d{i}        = [d{i}, tmp(ctmp)];
%     ddd         = [ddd, d{i}];
end
% x_est = x{1};
for i = 1:hr
    x{i} = unique(x{i});
    d{i} = unique(d{i});
%     x_est= [x_est; x{i}];
%     x_est= unique(x_est);
end
x_est = x{1};
for i = 2:hr
    x_est = union( x_est, x{i} );
end
x_unest = setdiff([1:n],x_est)';
size(A);
length(x_est);
length(x_unest);

%% Sub-Systems

F   = cell(hr,1);                 % Local A
Fl  = cell(hr,1);
H   = cell(hr,1);                 % Local C
D   = cell(hr,1);                 % Local Internal Input 
G   = cell(hr,1);                 % Local External Input
T   = cell(hr,1);
Q   = cell(hr,1);
R   = cell(hr,1);
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
% % %     H{i}    = C( i,x{i} );
    H{i}    = C( i,: );
    Q{i}    = QQ( x{i},x{i} );
    R{i}    = RR( i,i );
end

%% Prediction

Slff = cell(hr,1);
Slfd = cell(hr,1);
Sldd = cell(hr,1);
for i = 1:hr
    Slff{i} = Sk1k1( x{i},x{i} );            % Sff(k-1/k-1)_l
    Slfd{i} = Sk1k1( x{i},d{i} );            % Sfd(k-1/k-1)_l
    Sldd{i} = Sk1k1( d{i},d{i} );            % Sdd(k-1/k-1)_l
end

% S(k/k-1)_l
Slkk1 = cell(hr,1);
Zlkk1 = cell(hr,1);
xlkk1 = cell(hr,1);
tic
parfor i = 1:hr
    xlkk1{i} = xkk1(x{i});
    Slkk1{i} = F{i}*Slff{i}*F{i}' + F{i}*Slfd{i}*D{i}' + (F{i}*Slfd{i}*D{i}')' + D{i}*Sldd{i}*D{i}' + Q{i};
end
toc

% Z(k/k-1)_l
tic
% Zlkk1 = bandinv(Slkk1,L,x,hr,hc,0);
toc

% z(k/k-1)_l
zlkk1 = cell(hr,1);
tic
Stmp = A*Sk1k1*A' + QQ;
Ztmp = inv(Stmp);
toc
parfor i = 1:hr
%     Zlkk1{i} = Ztmp(x{i},x{i});
    Zlkk1{i} = inv(Slkk1{i});
    zlkk1{i} = Zlkk1{i}*xlkk1{i};
end

% Zkk1        = Ztmp;
% zkk1        = Zkk1*xkk1;
% ik          = C'*inv(RR)*y;
% Ik          = C'*inv(RR)*C;
% zkk         = zkk1 + ik;
% Zkk         = Zkk1 + Ik;
% 
% Zlkk = cell(hr,1);
% zlkk = cell(hr,1);
% parfor i = 1:hr
%     Zlkk{i} = Zkk(x{i},x{i});
%     zlkk{i} = zkk(x{i});
% end

%% Information Matrix (Local)

inff = cell(hr,1);
INFF = cell(hr,1);
parfor i = 1:hr
    inff{i} = H{i}(:,x{i})'*inv(R{i})*y(i);
    INFF{i} = H{i}(:,x{i})'*inv(R{i})*H{i}(:,x{i});
end

%% Local Update

Zlkk = cell(hr,1);
zlkk = cell(hr,1);
parfor i = 1:hr
    Zlkk{i} = Zlkk1{i} + INFF{i};
    zlkk{i} = zlkk1{i} + inff{i};
end

%% Kalman Filter

% tic
% Slkk = bandinv(Zlkk,L,x,hr,hc,0);
% toc

% P_yy    = cell(hr,1);
% P_xy    = cell(hr,1);
% P_xy    = zeros(hc,hr);
% P_yy    = zeros(hr,hr);
% K       = zeros(hc,hr);
% xlkk    = cell(hr,1);
% dy      = cell(hr,1);
% Slkk    = cell(hr,1);
% for i = 1:hr
%     dy{i}           = y(i) - H{i}*xkk1;
%     P_yy(i,i)       = RR(i,i) + H{i}(:,x{i})*Slkk1{i}*H{i}(:,x{i})';
%     P_xy(x{i},i)    = Slkk1{i}*H{i}(:,x{i})';
%     K(x{i},i)       = P_xy(x{i},i)*pinv(P_yy(i,i));
%     xlkk{i}         = xlkk1{i} + K(x{i},i)*dy{i};
%     Slkk{i}         = ( eye( size( K(x{i},i)*H{i}(:,x{i}) ) ) - K(x{i},i)*H{i}(:,x{i}) )*Slkk1{i};
% end

%% Data Fusion

% CI

% z = cell(hr,1);
% Z = cell(hr,1);
% for i = 1:hr
%     z{i} = zeros(n,1);
%     z{i}(x{i}) = zlkk{i};
%     Z{i} = zeros(n,n);
%     Z{i}(x{i},x{i}) = Zlkk{i};
% end
% zz = zeros(n,1);
% ZZ = zeros(n,n);
% for i = 1:hr
%     zz = zz + z{i};
%     ZZ = ZZ + Z{i};
% end
% Pkk             = pinv(ZZ);
% xkk             = Pkk*zz;
% xkk(x_unest)    = xkk1(x_unest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CI,EI,ICI - Information Matrix

tic
nnn = hr;
zf = zlkk;
Pf = Zlkk;
xf = x;
if type ~= 4 
    type;
    zfkk = cell(ceil(hr/2),1);
    Pfkk = cell(ceil(hr/2),1);
    xfkk = cell(ceil(hr/2),1);

    if nnn ~= 1
        if rem(nnn,2) == 0
            nnn = nnn;
            flag = 0;
        else
            nnn = nnn - 1;
            flag = 1;
        end
        nnn = nnn/2;
        while nnn >= 1
            ii = 0;
            parfor i = 1:(nnn)
                ii = i + (i - 1);
                [zfkk{i}, Pfkk{i}, xfkk{i}]  = fuze2(zf{ii},zf{ii+1},Pf{ii},Pf{ii+1},xf{ii},xf{ii+1},type);
            end
            if flag == 1
                loc = (nnn);
                [zfkk{loc}, Pfkk{loc}, xfkk{loc}]  = ...
                        fuze2(zfkk{loc},zf{2*nnn+1},Pfkk{loc},Pf{2*nnn+1},xfkk{loc},xf{2*nnn+1},type);
            end
            if nnn == 1
                nnn     = nnn/2;
                zkk     = zfkk{1};
                Pkk     = Pfkk{1};
                xfkk    = xfkk{1};  
            else
                zf = zfkk;
                Pf = Pfkk;
                xf = xfkk;
                zfkk = cell(ceil(nnn/2),1);
                Pfkk = cell(ceil(nnn/2),1);
                xfkk = cell(ceil(nnn/2),1);
                if ((rem(nnn/2,2) == 0)) || (nnn == 2)
                    nnn = nnn/2;
                    flag = 0;
                else
                    nnn = ((nnn) - 1)/2;
                    flag = 1;
                end
            end
            nnn = nnn;
        end
    end
    toc
    Ptmp            = pinv(Pkk);
    xtmp            = Ptmp*zkk;
else
    4;
    [xtmp,Ptmp] = fuze(zf,Pf,xf,hr,n,x_est);
end

Pkk             = zeros(n,n);
Pkk(x_est,x_est)= Ptmp;
xkk             = zeros(n,1);
xkk(x_est)      = xtmp;
xkk(x_unest)    = xkk1(x_unest);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EI - Co-Variance

% tic
% nnn = hr;
% % zf = zlkk;
% zf = cell(hr,1); 
% for i = 1:hr
%     zf{i} = Slkk{i}*zlkk{i};
% end
% Pf = Slkk;
% xf = x;
% zfkk = cell(ceil(hr/2),1);
% Pfkk = cell(ceil(hr/2),1);
% xfkk = cell(ceil(hr/2),1);
% 
% if nnn ~= 1
%     if rem(nnn,2) == 0
%         nnn = nnn;
%         flag = 0;
%     else
%         nnn = nnn - 1;
%         flag = 1;
%     end
%     nnn = nnn/2;
%     while nnn >= 1
%         ii = 0;
%         for i = 1:(nnn)
%             ii = i + (i - 1); 
% %             [X, X_a, X_b, Xij, xfkk{i}, xab_size, T1, T2, ia1, ia2, cx] ...
% %                         = fuse_cov(Pf{ii},Pf{ii+1},xf{ii},xf{ii+1});
% %             [zfkk{i}, Pfkk{i}] = fuse_mean(zf{ii},zf{ii+1}, X, X_a, X_b, Xij, xfkk{i}, xab_size, T1, T2, ia1, ia2, cx);
%             [zfkk{i}, Pfkk{i}, xfkk{i}]  = fuse2(zf{ii},zf{ii+1},Pf{ii},Pf{ii+1},xf{ii},xf{ii+1});
%         end
%         if flag == 1
%             loc = (nnn);
% %             [X, X_a, X_b, Xij, xfkk{loc}, xab_size, T1, T2, ia1, ia2, cx] ...
% %                         = fuse_cov(Pfkk{loc},Pf{2*nnn+1},xfkk{loc},xf{2*nnn+1});
% %             [zfkk{i}, Pfkk{i}] = fuse_mean(zfkk{loc},zf{2*nnn+1}, X, X_a, X_b, Xij, xfkk{loc}, xab_size, T1, T2, ia1, ia2, cx);
%             [zfkk{loc}, Pfkk{loc}, xfkk{loc}]  = ...
%                     fuse2(zfkk{loc},zf{2*nnn+1},Pfkk{loc},Pf{2*nnn+1},xfkk{loc},xf{2*nnn+1});
%         end
%         if nnn == 1
%             nnn     = nnn/2;
%             zkk     = zfkk{1};
%             Pkk     = Pfkk{1};
%             xfkk    = xfkk{1};  
%         else
%             zf = zfkk;
%             Pf = Pfkk;
%             xf = xfkk;
%             zfkk = cell(ceil(nnn/2),1);
%             Pfkk = cell(ceil(nnn/2),1);
%             xfkk = cell(ceil(nnn/2),1);
%             if ((rem(nnn/2,2) == 0)) || (nnn == 2)
%                 nnn = nnn/2;
%                 flag = 0;
%             else
%                 nnn = ((nnn) - 1)/2;
%                 flag = 1;
%             end
%         end
%         nnn = nnn;
%     end
% end
% toc
% % xkk = Pkk*zkk;
% xkk = zkk;
% Pkk = Pkk;
% size(xkk);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Skk1(x_est,x_est) = Pkk;
% xx = zeros(n,1);
% xx(x_est) = xkk;
% xx(x_unest) = xkk1(x_unest);
% 
% Pkk = Skk1;
% xkk = xx;

% if nnn ~= 1
%     if rem(nnn,2) == 0
%         nnn = nnn;
%         flag = 0;
%     else
%         nnn = nnn - 1;
%         flag = 1;
%     end
%     while nnn ~= 1
%         ii = 0;
%         parfor i = 1:(nnn/2)
%             if rem(i,2) == 0
%                 ii = i + 1;
%             else
%                 ii = i;
%             end
%             [zfkk{i}, Pfkk{i}, xfkk{i}]  = fuse2(zf{ii},zf{ii+1},Pf{ii},Pf{ii+1},xf{ii},xf{ii+1});
%             size(xfkk{i})
%         end
%         if flag == 1
%             loc = (nnn/2);
%             [zfkk{loc}, Pfkk{loc}, xfkk{loc}]  = ...
%                     fuse2(zfkk{loc},zf{nnn+1},Pfkk{loc},Pf{nnn+1},xfkk{loc},xf{nnn+1});
%         end
%         if nnn == 2
%             nnn     = nnn/2;
%             zkk     = zfkk{1};
%             Pkk     = Pfkk{1};
%             xfkk     = xfkk{1};  
%         else
%             zf = zfkk;
%             Pf = Pfkk;
%             xf = xfkk;
%             zfkk = cell(ceil(nnn/2),1);
%             Pfkk = cell(ceil(nnn/2),1);
%             xfkk = cell(ceil(nnn/2),1);
%             if (rem(nnn/2,2) == 0)
%                 nnn = nnn/2;
%                 flag = 0;
%             else
%                 nnn = ((nnn/2) - 1)/2;
%                 flag = 1;
%             end
%         end
%     end
% end
% toc
% xkk = Pkk*zkk;

% [cc,ia]                 = setdiff( [1:n],xfkk );
% T_tmp                   = eye(n,n);
% T_tmp(cc,cc)            = zeros(length(cc),length(cc));
% T_tmp                   = T_tmp([1:n],xfkk);
% T_tmp                   = sparse(T_tmp);
% xkk                     = T_tmp*xkk;
% Pkk                     = T_tmp*Pkk*T_tmp';





% % % % tic
% % % % xf      = x{1};
% % % % T_tmp   = cell(hr,1);
% % % % zf      = zlkk{1};
% % % % zfkk    = zlkk;
% % % % X      = Slkk{1};
% % % % % if rem(hr,2) == 0
% % % % for i = 1:hr-1
% % % %     x_tmp               = union( xf,x{i+1} );
% % % %     [cc,ia]             = setdiff( x_tmp,xf );
% % % %     T_tmp{i}            = eye(length(x_tmp));
% % % %     T_tmp{i}(ia,ia)     = zeros(length(ia),length(ia));
% % % %     zf                  = T_tmp{i}*zf;
% % % %         
% % % %     clear cc ia
% % % %     [cc,ia]             = setdiff( x_tmp,x{i+1} );
% % % %     T_tmp{i+1}          = eye(length(x_tmp));
% % % %     T_tmp{i+1}(ia,ia)   = zeros(length(ia),length(ia));
% % % %     zfkk{i+1}           = T_tmp{i+1}*zlkk{i+1};
% % % %         
% % % %     clear cc ia
% % % %     ia                  = find(~zf);
% % % %     ib                  = find(~zfkk{i+1});
% % % %     zf(ia)              = zfkk{i+1}(ia);
% % % %     zfkk{i+1}(ib)       = zf(ib);
% % % %     
% % % %     [Sa, Da] = eig(T_tmp{i}*X*T_tmp{i}');
% % % %     [Sb, Db] = eig(T_tmp{i+1}*Slkk{i+1}*T_tmp{i+1}');
% % % % 
% % % %     for i = 1:length(x_tmp)
% % % %         for j =1:3
% % % %             if Db(i,j) == 0
% % % %                 Gb(i,j) = 1;
% % % %             else
% % % %                 Gb(i,j) = Db(i,j);
% % % %             end
% % % %         end
% % % %     end
% % % % 
% % % %     for i = 1:3
% % % %         for j =1:3
% % % %             if Da(i,j) == 0
% % % %                 Ga(i,j) = 1;
% % % %             else
% % % %                 Ga(i,j) = Da(i,j);
% % % %             end
% % % %         end
% % % %     end
% % % % 
% % % %     Xb = Sa*(Ga^0.5)*Sb*Gb*inv(Sb)*(Ga^0.5)*inv(Sa);
% % % %     Xa = Sb*(Gb^0.5)*Sa*Ga*inv(Sa)*(Gb^0.5)*inv(Sb);
% % % % 
% % % %     X = inv(inv(Xa) + inv(Xb))
% % % %     zf = X*(inv(Xa)*zf + inv(Xb)*zfkk{i+1})
% % % %     xf = x_tmp;
% % % %     clear ia ib
% % % % end
% % % % Pkk = X;
% % % % xkk = Pkk*zf;
% % % % toc

%% Information Matrix (Union)

% infl = inff;
% INFl = INFF;
% for i = 1:hr
%     for j = [1:i-1 i+1:hr]
%         for z = intersect(x{i},x{j})
%             iz = find(x{i}==z);
%             jz = find(x{j}==z);
%             infl{i}(iz) = inff{i}(iz) + inff{j}(jz);
%         end
%     end
% end
% 
% for i = 1:hr
%     for j = [1:i-1 i+1:hr]
%         iz = []; jz = [];
%         kz = 1;
%         for z = intersect(x{i},x{j})
%             iz(kz) = find(x{i}==z);
%             jz(kz) = find(x{j}==z);
%             kz = kz + 1;
%         end
%         INFl{i}(iz,iz) = INFF{i}(iz,iz) + INFF{j}(jz,jz);
% %         end
%     end
% end

%% Update Equation

% % z(k/k)_l & Z(k/k)_l
% zlkk = cell{hr,1};
% Zlkk = cell{hr,1};
% for i = 1:hr
%     zlkk{i} = zlkk1{i} + infl{i};
%     Zlkk{i} = Zlkk1{i} + INFl{i};
% end
% 
% % z(k/k) & Z(k/k)
% for i = 1:hr
%     zkk( x{i} ) = zlkk{i};
%     Zkk( x{i},x{i} ) = Zlkk{i};       
% end
% 
% % x(k/k)
% Skk = inv(Zkk);
% [bandl,bandu]=bandwidth(Zkk);
% % Skk = bandinv(Zkk,max(bandl,bandu));
% xkk = Skk*zkk;




% % clear iT dx loc1 loc2 iii jjj crossmat_temp iturb
% %     
% % % Secondly, calculate the autocorrelation of output
% % rho_locl.auto = sparse(size(outputLocArray,1),size(outputLocArray,1));
% % for(iii = 1:size(outputLocArray,1)) % for each output
% %     loc1 = outputLocArray(iii,:);
% %     for jjj = iii:strucObs.nrobs % for each output
% %         loc2 = outputLocArray(jjj,:);
% %         dx = sqrt(sum((loc1-loc2).^2));
% %         rho_locl.auto(iii,jjj) = localizationGain( dx, strucObs.f_locl, strucObs.l_locl );
% %     end;
% % end;








