function [xkk Pkk Ptmp Slkk K] = distributed_linear2( strucObs,x,d,p,hr,n, F,D,E,G,H,Q,R, yy,yyy, xkk1,Slkk, x_est,x_unest, P_unest, type,typeCZ, weight,constant,K  );
% [xkk Pkk] = distributed_linear( x,d,F,D,G,H,Q,R, y, xkk1,xk1k1,Sk1k1, type );


% S(k/k-1)_l
xlkk1 = cell(hr,1);
for i = 1:hr
    xlkk1{i}    = xkk1(x{i});

    y{i}        = yy( : );
    
end


%% Kalman Filter (Update)

if strcmp(typeCZ,'C')
    xlkk    = cell(hr,1);
    dy      = cell(hr,1);
    parfor i = 1:hr
        dy{i}           = y{i} - H{i}*xlkk1{i};
        xlkk{i}         = xlkk1{i} + K{i}*dy{i};
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
if strcmp(type,'CIN')||strcmp(type,'CI2')||strcmp(type,'EI')||strcmp(type,'ICI')
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
            for i = 1:(nnn)
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
elseif strcmp(type,'CI')
    for i = 1:hr
        tmp             = zeros(n,n);
        tmp(x{i},x{i})  = Pf{i};
        Pf{i}           = tmp(x_est,x_est);
        tmp             = zeros(n,1);
        tmp(x{i})       = zf{i};
        zf{i}           = tmp(x_est);
    end  
    A = ones(1,hr);B = 1;
    omega0 = (1/hr)*ones(1,hr);
    options = optimset('tolx',1.0e-8,'MaxFunEvals',25,'Display','off');
    omega = fmincon(@f,omega0',A,B,[],[],0*A,A,[],options,Pf,xf,hr);
    %omega = omega0;
    l_est = length(x_est);
    Ptmp = zeros(l_est,l_est);
    xtmp = zeros(l_est,1);
    for i = 1:hr
        Ptmp = Ptmp + omega(i)*Pf{i};
        xtmp = xtmp + omega(i)*zf{i};
    end
    Ptmp = pinv(Ptmp);
    xtmp = Ptmp*xtmp;
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

