function S = bandinv(Z,L,x,hr,hc,type);
% S = bandinv(Z,L,x,hr,hc,type);
%   type   : 1 if Z is a matrix
%            0 if Z is a cell
%   Z      : Matrix to be Inverted
%   L      : Bandwidth
%   x      : Local States
%   hr     : Number of Subsystems
%   hc     : Total number of states in a Large-Scale System

% L       = 1;
% hc      = 5;
% Z       = bandmatrix(hc,L,L);
% hr      = 3;
% x       = cell(hr,1);
% x{1}    = [1:3];
% x{2}    = [2:4];
% x{3}    = [4:5];
% S = bandinv(Z,L,x,hr,hc,type);

% clear all
% close all
% clc
% 
% type    = 1;
% L       = 3;
% hc      = 5;
% hr      = 3; 
% Z       = bandmatrix(hc,L,L);
% ZZ      = cell(3,1);
% T       = cell(3,1);
% x       = cell(3,1);
% x{1}    = [1:3];
% x{2}    = [2:4];
% x{3}    = [4:5];
% for i = 1:3
%     T{i}            = zeros(n,n);
%     l               = length(x{i});
%     T{i}(x{i},x{i}) = eye(l,l);
% end
% for i = 1:3
% %     ZZ{i}   = Z(x{i},x{i});
%     ZZ{i} = T{i}*Z*T{i}';
% end

ZZ = cell(hr,1);
for i = 1:hr
%     ZZ{i} = zeros(hc,hc);
    if type == 1
%         ZZ{i}(x{i},x{i}) = Z(x{i},x{i});
        ZZ{i} = Z(x{i},x{i});
    else
%         ZZ{i}(x{i},x{i}) = Z{i};
        ZZ{i} = Z{i};
    end
end

M = cell(hr,1);
m = cell(hr,1);
P = cell(hr,1);
S = cell(hr,1);
tic
parfor i=1:hr
    M{i} = diag(diag(ZZ{i}));
%     m{i} = zeros(hc,hc);
%     for ii = 1:hc
%         if M{i}(ii,ii) ~= 0
%             m{i}(ii,ii) = inv(M{i}(ii,ii));
            m{i} = diag(1./diag(ZZ{i}));
%         end
%     end
    P{i} = m{i}*(ZZ{i} - M{i});
end
toc
% inv(diag(diag(Z)))*(Z - diag(diag(Z)));
[rZ cZ] = size(Z);

% for i = 1:hr
%     M{i} = sparse(M{i});
%     m{i} = sparse(m{i});
%     P{i} = sparse(P{i});
% %     ZZ{i}= sparse(ZZ{i});
%     S{i} = sparse(S{i});
% end

%% Working

% parfor ij = 1:hr
%     S{ij}       = inv(ZZ{ij});
%     tmp         = S{ij};
%     k           = 1;
%     while ((k == 1) || (abs( norm(eig(tmp)-eig(S{ij}),Inf) ) >= 10e-10)) && k < 20
%         tmp = S{ij};
%         ik = 1;
%         l = length(x{ij});
%         for i = x{ij}
%             jk = 1;
%             for j = x{ij}
%                 if ((abs(i - j) <= L) | (jk == 1)) | (ik == max(x{ij}))
%                     S{ij}(ik,jk) = - P{ij}(ik,:)*S{ij}(:,jk) + m{ij}(ik,jk);
%                 else
%                     S{ij}(ik,jk) = S{ij}(ik,jk-1)*pinv(S{ij}(ik+1,jk-1))*S{ij}(ik+1,jk);
%                 end
%                 jk = jk + 1;
%             end
%             ik = ik + 1;
%         end
%         k = k + 1;
%     end
% end

%% L-Banded Inverse [ZS = b] (Not Working)
% tic
% for ij = 1:hr
%     S{ij}       = m{ij};
% %     S{ij}       = inv(ZZ{ij});
%     tmp         = S{ij};
%     k           = 1;
%     while ((k == 1) || (abs( norm(eig(tmp)-eig(S{ij}),Inf) ) >= 10e-10)) && k < 20
%         tmp = S{ij};
%         ik = 1;
%         l = length(x{ij});
%         S{ij} = - P{ij}*S{ij} + m{ij};
%         for i = 1:l
%             i = x{ij}(i);
%             jk = 1;
%             for j = 1:l
%                 j = x{ij}(j);
%                 if ((abs(i - j) > L) & (jk ~= 1)) & (ik ~= l)
% %                     S{ij}(ik,jk) = S{ij}(ik,jk-1)*pinv(S{ij}(ik+1,jk-1))*S{ij}(ik+1,jk);
%                     S{ij}(ik,jk) = - P{ij}(ik,:)*S{ij}(:,jk) + m{ij}(ik,jk);
%                 end
%                 jk = jk + 1;
%             end
%             ik = ik + 1;
%         end
%         k = k + 1;
%     end
% end
% toc
%% De-Centralized (Substates inorder)

parfor ij = 1:hr
    S{ij}       = m{ij};
%     S{ij}       = inv(ZZ{ij});
    tmp         = S{ij};
    k           = 1;
    if ij == 1
        xtmp = union(x{ij},x{ij+1});
        l = max(xtmp);
        Stmp = zeros(l,l);
        Ptmp = zeros(l,l);
        Stmp(x{ij},x{ij}) = m{ij};
        Stmp(x{ij+1},x{ij+1}) = m{ij+1};
        Ptmp(x{ij},x{ij}) = P{ij};
        Ptmp(x{ij+1},x{ij+1}) = P{ij+1};
        mtmp = Stmp;
    elseif (ij>1)&&(ij<hr)
        xtmp = union(x{ij},x{ij-1});
        xtmp = union(xtmp,x{ij+1});
        l = max(xtmp);
        Stmp = zeros(l,l);
        Ptmp = zeros(l,l);
        Stmp(x{ij-1},x{ij-1}) = m{ij-1};
        Stmp(x{ij},x{ij}) = m{ij};
        Stmp(x{ij+1},x{ij+1}) = m{ij+1};
        Ptmp(x{ij-1},x{ij-1}) = P{ij-1};
        Ptmp(x{ij},x{ij}) = P{ij};
        Ptmp(x{ij+1},x{ij+1}) = P{ij+1};
        mtmp = Stmp;
    else
        xtmp = union(x{ij},x{ij-1});
        l = max(xtmp);
        Stmp = zeros(l,l);
        Ptmp = zeros(l,l);
        Stmp(x{ij-1},x{ij-1}) = m{ij-1};
        Stmp(x{ij},x{ij}) = m{ij};
        Ptmp(x{ij-1},x{ij-1}) = P{ij-1};
        Ptmp(x{ij},x{ij}) = P{ij};
        mtmp = Stmp;
    end
    while ((k == 1) || (abs( norm(eig(tmp)-eig(Stmp),Inf) ) >= 10e-10)) && k < 20
        tmp = Stmp;
        ik = 1;
        for i = 1:l
            i = xtmp(ij);
            jk = 1;
            for j = 1:l
                j = xtmp(ij);
                if ((abs(i - j) <= L) | (jk == 1)) | (ik == j)
                    Stmp(ik,jk) = - Ptmp(ik,:)*Stmp(:,jk) + mtmp(ik,jk);
                else
                    Stmp(ik,jk) = Stmp(ik,jk-1)*pinv(Stmp(ik+1,jk-1))*Stmp(ik+1,jk);
                end
                jk = jk + 1;
            end
            ik = ik + 1;
        end
        k = k + 1;
    end
    S{ij} = Stmp(x{ij},x{ij});
end

%% Distributed

% Stmp = zeros(hc,hc);
% Ptmp = zeros(hc,hc);
% for ij = 1:hr
%     Stmp(x{ij},x{ij}) = m{ij};
%     Ptmp(x{ij},x{ij}) = P{ij};
% end
% mtmp    = Stmp;
% 
% %     Stmp = zeros(hc,hc);
% %     Ptmp = zeros(hc,hc);
% %     Stmp(x{ij},x{ij}) = m{ij};
% %     Ptmp(x{ij},x{ij}) = P{ij};
% %     mtmp    = Stmp;
% 
% for ij = 1:hr
%     k       = 1;
%     tmp     = Stmp;
%     while ((k == 1) || (abs( norm(eig(tmp)-eig(Stmp),Inf) ) >= 10e-10)) && k < 20
%         tmp = Stmp;
%         ik = 1;
%         l = length(x{ij});
%         for i = 1:l
%             i = x{ij}(i);
%             jk = 1;
%                 for j = 1:l
%                     i = x{ij}(j);
%                     if ((abs(i - j) <= L) | (jk == 1)) | (ik == hc)
%                         Stmp(ik,jk) = - Ptmp(ik,:)*Stmp(:,jk) + mtmp(ik,jk);
%                     else
%                         Stmp(ik,jk) = Stmp(ik,jk-1)*pinv(Stmp(ik+1,jk-1))*Stmp(ik+1,jk);
%                     end
%                     jk = jk + 1;
%                 end
%                 ik = ik + 1;
%             end
%             k = k + 1;
%     end
%     S{ij} = Stmp(x{ij},x{ij});
% end
