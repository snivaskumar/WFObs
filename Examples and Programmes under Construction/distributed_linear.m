function [xkk Pkk] = distributed_linear( x,d,p,hr,n, F,D,G,H,Q,R, yy, xkk1,xk1k1,Sk1k1, x_est,x_unest, P_unest, type,typeCZ );
% [xkk Pkk] = distributed_linear( x,d,F,D,G,H,Q,R, y, xkk1,xk1k1,Sk1k1, type );
% tic
% xkk1    = xkk1(p);
% xk1k1   = xk1k1(p);
% Sk1k1   = Sk1k1(p,p);
% [Y,I] = sort(p);
% toc

% Prediction
Slff = cell(hr,1);
Slfd = cell(hr,1);
Sldd = cell(hr,1);
% tic
for i = 1:hr
    Slff{i} = Sk1k1( x{i},x{i} );            % Sff(k-1/k-1)_l
    Slfd{i} = Sk1k1( x{i},d{i} );            % Sfd(k-1/k-1)_l
    Sldd{i} = Sk1k1( d{i},d{i} );            % Sdd(k-1/k-1)_l
end
% toc

% S(k/k-1)_l
Slkk1 = cell(hr,1);
Zlkk1 = cell(hr,1);
xlkk1 = cell(hr,1);
% tic
parfor i = 1:hr
    xlkk1{i}    = xkk1(x{i});
    Slkk1{i}    = F{i}*Slff{i}*F{i}'...
                + F{i}*Slfd{i}*D{i}'...
                + (F{i}*Slfd{i}*D{i}')'...
                + D{i}*Sldd{i}*D{i}' + Q{i};
    y{i}        = yy( : );
end
% toc

%% Information Filter (Update)

% Z(k/k-1)_l, z(k/k-1)_l
zlkk1 = cell(hr,1);
parfor i = 1:hr
    Zlkk1{i} = inv(Slkk1{i});
    zlkk1{i} = Zlkk1{i}*xlkk1{i};
end

if typeCZ == 2
    inff = cell(hr,1);
    INFF = cell(hr,1);
    parfor i = 1:hr
        inff{i} = H{i}'*inv(R{i})*y{i};
        INFF{i} = H{i}'*inv(R{i})*H{i};
    end
    
    Zlkk = cell(hr,1);
    zlkk = cell(hr,1);
    parfor i = 1:hr
        Zlkk{i} = Zlkk1{i} + INFF{i};
        zlkk{i} = zlkk1{i} + inff{i};
    end
end

%% Kalman Filter (Update)

if typeCZ == 1
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
                if typeCZ == 1
                    [zfkk{i}, Pfkk{i}, xfkk{i}]  = fuse2(zf{ii},zf{ii+1},Pf{ii},Pf{ii+1},xf{ii},xf{ii+1},type);
                else
                    [zfkk{i}, Pfkk{i}, xfkk{i}]  = fuze2(zf{ii},zf{ii+1},Pf{ii},Pf{ii+1},xf{ii},xf{ii+1},type);
                end
            end
            if flag == 1
                loc = (nnn);
                if typeCZ == 1
                    [zfkk{loc}, Pfkk{loc}, xfkk{loc}]  = ...
                        fuse2(zfkk{loc},zf{2*nnn+1},Pfkk{loc},Pf{2*nnn+1},xfkk{loc},xf{2*nnn+1},type);
                else
                    [zfkk{loc}, Pfkk{loc}, xfkk{loc}]  = ...
                            fuze2(zfkk{loc},zf{2*nnn+1},Pfkk{loc},Pf{2*nnn+1},xfkk{loc},xf{2*nnn+1},type);
                end
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
    if typeCZ == 1
        Ptmp            = Pkk;
        xtmp            = zkk;
    else
        Ptmp            = pinv(Pkk);
        xtmp            = Ptmp*zkk;
    end
else
    4;
    [xtmp,Ptmp] = fuze(zf,Pf,xf,hr,n,x_est,typeCZ);
end

Pkk             = 5*eye(n,n);
Pkk(x_est,x_est)= Ptmp;

xkk             = zeros(n,1);
xkk(x_est)      = xtmp;
xkk(x_unest)    = xkk1(x_unest);

% Pkk = Pkk(I,I);
% xkk = xkk(I);
% toc

