function [xkk Pkk omega] = distributed_linear( strucObs,x,d,p,hr,n, F,D,E,G,H,Q,R, yy,yyy, xkk1,xk1k1,Sk1k1, x_est,x_unest, P_unest, type,typeCZ, weight,constant );
% [xkk Pkk] = distributed_linear( x,d,F,D,G,H,Q,R, y, xkk1,xk1k1,Sk1k1, type );
% tic
% xkk1    = xkk1(p);
% xk1k1   = xk1k1(p);
% Sk1k1   = Sk1k1(p,p);
% [Y,I] = sort(p);
% toc
omega = [];
% Prediction
Slff = cell(hr,1);
Slfd = cell(hr,1);
Sldd = cell(hr,1);

% Slee = cell(hr,1);
% Slfe = cell(hr,1);
% Slde = cell(hr,1);
if strucObs.Optimize == 0
    Slee_tmp = Sk1k1( x_unest,x_unest );
else
    Slee_tmp = sparse(diag(diag( Sk1k1( x_unest,x_unest ) )));
end
% tic
for i = 1:hr
    Slff{i} = Sk1k1( x{i},x{i} );            % Sff(k-1/k-1)_l
    Slfd{i} = Sk1k1( x{i},d{i} );            % Sfd(k-1/k-1)_l
    Sldd{i} = Sk1k1( d{i},d{i} );            % Sdd(k-1/k-1)_l
    if strucObs.extremeOptimize == 1
        factor = strucObs.superOptimizeFactor;
        indices         = find( abs(Sldd{i})<factor );
        Sldd{i}(indices)= 0;
        Sldd{i}         = sparse(Sldd{i});
        
        indices         = find( abs(Slfd{i})<factor );
        Slfd{i}(indices)= 0;
        Slfd{i}         = sparse(Slfd{i});
        
        indices         = find( abs(Slff{i})<factor );
        Slff{i}(indices)= 0;
        Slff{i}         = sparse(Slff{i});
    end
    Slee{i} = Slee_tmp;
    if strucObs.Optimize == 0
%     Unoptimized (Entire matrices):
        Slfe{i} = Sk1k1( x{i},x_unest );
        Slde{i} = Sk1k1( d{i},x_unest );
    else        
%     Optimized (Only the diagonal entries):
        len_x       = length(x{i});
        len_d       = length(d{i});
        len_xunest  = length(x_unest);

%         Slee{i} = sparse(diag(diag(Slee_tmp)));
        Slfe{i} = sparse(zeros(len_x,len_xunest));
        Slfe{i}(1:min(len_x,len_xunest),1:min(len_x,len_xunest)) = sparse(diag(diag(Sk1k1( x{i},x_unest ))));
        Slde{i} = sparse(zeros(len_d,len_xunest));
        Slde{i}(1:min(len_d,len_xunest),1:min(len_d,len_xunest)) = sparse(diag(diag(Sk1k1( d{i},x_unest ))));
    end
end
% toc

typeOutput = upper( strucObs.typeOutput );
% S(k/k-1)_l
Slkk1 = cell(hr,1);
xlkk1 = cell(hr,1);
kk = cell(hr,1);
jj = cell(hr,1);
% tic
parfor i = 1:hr
    xlkk1{i}    = xkk1(x{i});
%     Slkk1{i}    = F{i}*Slff{i}*F{i}'...
%                 + F{i}*Slfd{i}*D{i}'...
%                 + (F{i}*Slfd{i}*D{i}')'...
%                 + D{i}*Sldd{i}*D{i}' + Q{i};
            
%     Slkk1{i}    = F{i}*Slff{i}*F{i}'...
%                 + F{i}*Slfd{i}*D{i}'...
%                 + (F{i}*Slfd{i}*D{i}')'...
%                 + D{i}*Sldd{i}*D{i}' + Q{i}...
%                 + F{i}*Slfe{i}*E{i}' + (F{i}*Slfe{i}*E{i}')'...
%                 + D{i}*Slde{i}*E{i}' + (D{i}*Slde{i}*E{i}')'...
%                 + E{i}*Slee{i}*E{i}';

    S1 = F{i}*Slff{i}*F{i}'; S2 = F{i}*Slfd{i}*D{i}';
    S3 = D{i}*Sldd{i}*D{i}';    
    S4 = D{i}*Slde{i}*E{i}'; S5 = E{i}*Slee{i}*E{i}'; 
    S6 = F{i}*Slfe{i}*E{i}';
    
%     Slkk1{i} = S1 + S2 + S2' + S3 + Q{i};    
    Slkk1{i} = S1 + S2 + S2'...
             + S3 + S4 + S4'...
             + S5 + S6 + S6' + Q{i};
            
    if strcmp(typeOutput,'ALL')  
        y{i}        = yy( : );
    else
        obs_array = strucObs.obs_array;
        len = length( obs_array );
        y{i} = [];
        kk{i} = [];
        jj{i} = [];
        for j = 1:len
            k = obs_array(j);
            if sum( x{i}==k )
                y{i} = [y{i}; yyy(k)];
                kk{i} = [kk{i};k];
                jj{i} = [jj{i};j];
            end
        end
        H{i} = H{i}(jj{i},:);
        [r,c] = size( H{i} );
        R{i} = strucObs.R_k*eye(r,r);
    end
    
end
% toc

%% Information Filter (Update)

if strcmp(typeCZ,'Z')
    Zlkk1 = cell(hr,1);     inff = cell(hr,1);
    zlkk1 = cell(hr,1);     INFF = cell(hr,1);
    Zlkk = cell(hr,1);
    zlkk = cell(hr,1);
    parfor i = 1:hr
        Zlkk1{i} = inv(Slkk1{i});
        zlkk1{i} = Zlkk1{i}*xlkk1{i};
        
        INFF{i} = H{i}'*inv(R{i})*H{i};
        inff{i} = H{i}'*inv(R{i})*y{i};
        
        Zlkk{i} = Zlkk1{i} + INFF{i};
        zlkk{i} = zlkk1{i} + inff{i};
    end
end

%% Kalman Filter (Update)

if strcmp(typeCZ,'C')
    Zlkk1{i} = Slkk1{i};
    P_yy    = cell(hr,1);
    P_xy    = cell(hr,1);
    K       = cell(hr,1);
    xlkk    = cell(hr,1);
    dy      = cell(hr,1);
    Slkk    = cell(hr,1);
    parfor i = 1:hr
        dy{i}           = y{i} - H{i}*xlkk1{i};
        P_yy{i}         = R{i} + H{i}*Slkk1{i}*H{i}';
        P_xy{i}         = Slkk1{i}*H{i}';
        K{i}            = P_xy{i}*pinv(P_yy{i});
        xlkk{i}         = xlkk1{i} + K{i}*dy{i};
        Slkk{i}         = ( eye( size( K{i}*H{i} ) ) - K{i}*H{i} )*Slkk1{i};
    end
    Zlkk = Slkk;
    zlkk = xlkk;
end

%% Data Fusion
% CI,EI,ICI,IFAC - Information Matrix

% tic
nnn = hr;
zf = zlkk;
Pf = Zlkk;
xf = x;
filter          = strucObs.filtertype; 
iteration       = strucObs.fusion_CIiteration; 
if strcmp(type,'CIN')||strcmp(type,'CI2')||strcmp(type,'CI')||strcmp(type,'EI')||strcmp(type,'ICI')
    zfkk = cell(ceil(hr/2),1);
    Pfkk = cell(ceil(hr/2),1);
    xfkk = cell(ceil(hr/2),1);
    omega = cell(ceil(hr/2),1);
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
                if strcmp(typeCZ,'C')
                    [zfkk{i}, Pfkk{i}, xfkk{i}, omega{i}]  = fuse2(zf{ii},zf{ii+1},Pf{ii},Pf{ii+1},xf{ii},xf{ii+1},type, weight,constant,iteration,filter);
                elseif strcmp(typeCZ,'Z')
                    [zfkk{i}, Pfkk{i}, xfkk{i}, omega{i}]  = fuze2(zf{ii},zf{ii+1},Pf{ii},Pf{ii+1},xf{ii},xf{ii+1},type, weight,constant,iteration);
                end
            end
            if flag == 1
                loc = (nnn);
                if strcmp(typeCZ,'C')
                    [zfkk{loc}, Pfkk{loc}, xfkk{loc}, omega{loc}]  = ...
                        fuse2(zfkk{loc},zf{2*nnn+1},Pfkk{loc},Pf{2*nnn+1},xfkk{loc},xf{2*nnn+1},type, weight,constant,iteration,filter);
                elseif strcmp(typeCZ,'Z')
                    [zfkk{loc}, Pfkk{loc}, xfkk{loc}, omega{loc}]  = ...
                            fuze2(zfkk{loc},zf{2*nnn+1},Pfkk{loc},Pf{2*nnn+1},xfkk{loc},xf{2*nnn+1},type, weight,constant,iteration);
                end
            end
            if nnn == 1
                nnn     = nnn/2;
                zkk     = zfkk{1};
                Pkk     = Pfkk{1};
                xfkk    = xfkk{1};  
                omega   = omega{1};
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
    if strcmp(typeCZ,'C')
        Ptmp            = Pkk;
        xtmp            = zkk;
    elseif strcmp(typeCZ,'Z')
        Ptmp            = pinv(Pkk);
        xtmp            = Ptmp*zkk;
    end
elseif strcmp(type,'IFAC')
    4;
    typeIFAC    = strucObs.IFAC_type;
    typeWeight  = upper(strucObs.IFACWeight);
    weight      = strucObs.IFACConstantWeight;
    [xtmp,Ptmp] = fuze(zf,Pf,xf,hr,n,Zlkk1,xkk1,x_est,typeCZ,typeIFAC,typeWeight,weight);
% elseif strcmp(type,'CI')
%     for i = 1:hr
%         tmp             = zeros(n,n);
%         tmp(x{i},x{i})  = Pf{i};
%         Pf{i}           = tmp(x_est,x_est);
%         tmp             = zeros(n,1);
%         tmp(x{i})       = zf{i};
%         zf{i}           = tmp(x_est);
%     end  
%     A = ones(1,hr);B = 1;
%     omega0 = (1/hr)*ones(1,hr);
%     options = optimset('tolx',1.0e-8,'MaxFunEvals',25,'Display','off');
%     omega = fmincon(@f,omega0',A,B,[],[],0*A,A,[],options,Pf,xf,hr);
%     %omega = omega0;
%     l_est = length(x_est);
%     Ptmp = zeros(l_est,l_est);
%     xtmp = zeros(l_est,1);
%     for i = 1:hr
%         Ptmp = Ptmp + omega(i)*Pf{i};
%         xtmp = xtmp + omega(i)*zf{i};
%     end
%     Ptmp = pinv(Ptmp);
%     xtmp = Ptmp*xtmp;
end

Punest = strucObs.Punest;
if strcmp(type,'NO FUSION')||strcmp(type,'NO')     %No Fusion
    Pkk             = Punest*eye(n,n);
    xkk             = xkk1;
%     xkk             = zeros(n,1);
    for  i = 1:hr
        if strcmp(typeCZ,'C')
            xkk(xf{i}) = zf{i};
            Pkk(xf{i},xf{i}) = Pf{i};
        elseif strcmp(typeCZ,'Z')
            PPtmp = inv(Pf{i});
            Pkk(xf{i},xf{i}) = PPtmp;
            xkk(xf{i}) = PPtmp*zf{i};
        end
    end
%     xkk(x_unest)    = xkk1(x_unest);
else        %Fusion
    Pkk             = Punest*eye(n,n);
    Pkk(x_est,x_est)= Ptmp;

    xkk             = zeros(n,1);
    xkk(x_est)      = xtmp;
    xkk(x_unest)    = xkk1(x_unest);
end

% l_est = length( x_est );
% l_unest = length( x_unest );
% 
% Pest = Ptmp;
% Punest = strucObs.Punest*eye(l_unest,l_unest);
% 
% xest = xtmp;
% xunest = xkk1(x_unest);
% 
% [xkk, Pkk, ~, ~] = fuse2(xest,xunest,Pest,Punest,x_est,x_unest,'CIN',weight,constant,iteration,filter);

% toc

