clear all
close all
clc

hr = 4;
x       = cell(hr,1);
x{1}    = [1,2,3]';
x{2}    = [1,4]';
x{3}    = [4,5]';
x{4}    = [2,3,6]';

% zz = [1.8095; -0.1157; 2.8389; 1.1741; 1.7378; 1.8185];
zk1k1= [1.8121; -0.1938; 2.7738; 1.2965; 1.7856; 1.7763];
% load('/Users/Nivas_Kumar/Desktop/haha.mat','Zt1t1','ztt');
% Zlkk = Zt1t1;
% zlkk = ztt;
load('/Users/Nivas_Kumar/Desktop/haha1.mat','Ptmp','ztt');
Slkk = Ptmp;
zlkk = ztt;

% zlkk        = cell(hr,1);
% zlkk{1}     = [2.5914; -0.3437; 2.2813];
% zlkk{2}     = [0.7771; 2.1951];
% zlkk{3}     = [1.3846; 1.5312];
% zlkk{4}     = [0.7526; 1.9811; 1.8081];
% 
% Slkk        = cell(hr,1);
% Slkk{1}     = [0.6283   -0.3062   -0.3185;
%                -0.3062    0.6265   -0.3168;
%                -0.3185   -0.3168    0.6383];
% Slkk{2}     = [0.5081   -0.5031;
%                -0.5031    0.5082];
% Slkk{3}     = [0.4502   -0.4442;
%                -0.4442    0.4482];
% Slkk{4}     = [0.6240   -0.3378   -0.2825;
%                -0.3378    0.6026   -0.2613;
%                -0.2825   -0.2613    0.5468];

% hr = 5;
% % zz = [-1 10 -5.3 2.5 7.5];
% zlkk        = cell(hr,1);
% zlkk{1}     = [-1.13223 10.23131]';
% zlkk{2}     = [9.96432 -5.325672]';
% zlkk{3}     = [-5.324243 2.546464]';
% zlkk{4}     = [-0.989435 2.244234]';
% zlkk{5}     = [-5.024243 7.546464]';
% Slkk        = cell(hr,1);
% lx          = zeros(hr,1);
% for i = 1:hr
%     lx(i)       = length(zlkk{i});
%     Slkk{i}     = bandmatrix(lx(i),lx(i),lx(i));
% end
% 
% x       = cell(hr,1);
% x{1}    = [1,2]';
% x{2}    = [2,3]';
% x{3}    = [3,4]';
% x{4}    = [1,4]';
% x{5}    = [3,5]';

%%

tic
nnn = hr;
zf = zlkk;
Pf = Slkk;
xf = x;
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
        for i = 1:(nnn)
            ii = i + (i - 1); 
%             [X, X_a, X_b, Xij, xfkk{i}, xab_size, T1, T2, ia1, ia2, cx] ...
%                 = fuse_cov(Pf{ii},Pf{ii+1},xf{ii},xf{ii+1});
%             [zfkk{i}, Pfkk{i}] = fuse_mean(zf{ii},zf{ii+1}, X, X_a, X_b, Xij, xfkk{i}, xab_size, T1, T2, ia1, ia2, cx);

            [zfkk{i}, Pfkk{i}, xfkk{i}]  = fuse2(zf{ii},zf{ii+1},Pf{ii},Pf{ii+1},xf{ii},xf{ii+1});
%         zfkk{i}
        end
        if flag == 1
            loc = (nnn);
%             [X, X_a, X_b, Xij, xfkk{loc}, xab_size, T1, T2, ia1, ia2, cx] ...
%                 = fuse_cov(Pfkk{loc},Pf{2*nnn+1},xfkk{loc},xf{2*nnn+1});
%             [zfkk{loc}, Pfkk{loc}] = fuse_mean(zfkk{loc},zf{2*nnn+1}, X, X_a, X_b, Xij, xfkk{loc}, xab_size, T1, T2, ia1, ia2, cx);
            
            [zfkk{loc}, Pfkk{loc}, xfkk{loc}]  = ...
                    fuse2(zfkk{loc},zf{2*nnn+1},Pfkk{loc},Pf{2*nnn+1},xfkk{loc},xf{2*nnn+1});
%         zfkk{loc}
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
zkk
for i = 1:length(zkk)
    if zkk(i)<0
        haha(i) = -abs(zkk(i));
    else
        haha(i) = abs(zkk(i));
    end
end
haha = haha'
% xkk = Pkk*zkk;

%%

% H       = [];
% y       = [];
% Cv      = []; 
% kk      = 1; 
% w       = rand(hr,1);
% w       = w/sum(w);
% % w = (1/hr)*ones(hr,1);
% xx      = []; 
% for i = 1:hr
%     tmp     = double([1:6]==x{i});
%     H       = [H; tmp];
%     y       = [y; zlkk{i}];
%     ll      = length(Slkk{i});
%     Cv(kk:kk+ll-1,kk:kk+ll-1)= [eig(Slkk{i}).*eye(size(Slkk{i}))];
% %     ttmp(i) = sum(det(Slkk{i}));
% %     ttmp(i) = 1/ttmp(i);
%     ttmp(kk:kk+ll-1) = 1./(eig(Slkk{i}));
%     kk      = kk + ll;
%     xx      = [xx, x{i}']; 
% end
% xtmp = [];
% kk = 1;
% for i = 1:hr
%     xtmp(i) = sum(double(xx == i));
% %     [cc,ia] = find(xx,i)
%     ww = xtmp(i)/length(xx);
%     loc = find(double(xx == i).*ttmp);
%     w(loc) = ww.*(ttmp(loc)/sum(ttmp(loc)));
% end
% kk = 1;
% Cv = (1./w).*Cv;
% % for i = 1:hr
% %     Cv(kk:kk+ll-1,kk:kk+ll-1)   = (1/w(i))*Cv(kk:kk+ll-1,kk:kk+ll-1);
% %     kk      = kk + ll;
% % end
% xe = inv(H'*Cv*H)*H'*inv(Cv)*y
% 100*xe