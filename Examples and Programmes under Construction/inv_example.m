clear all
close all
clc
x = cell(4,1);
x{1}    = [1,2,3]';
x{2}    = [1,4]';
x{3}    = [4,5]';
x{4}    = [2,3,6]';
n = 6;
L = 1;
P = bandmatrix(n,L,L);
% P = rand(n,n);
P = (P*P')/norm(P*P');
% Z  = inv(P);
% P1 = P(x{1},x{1});
% Z1 = inv(P1);
% P2 = P(x{2},x{2});
% Z2 = inv(P2);
% P3 = P(x{3},x{3});
% Z3 = inv(P3);
% P4 = P(x{4},x{4});
% Z4 = inv(P4);

% PP = [P(x{1},:); P([4:6],x{1}),zeros(3,3)];

D = diag(diag(P));
Z = pinv(D);
ZZ = cell(100,1);
ztmp1 = inv(D);
ztmp2 = inv(D);
k = 1;
while (k==1)||( (norm( eig(ztmp1) - eig(ztmp2),Inf )>1e-10)&&(k<50) )
    ztmp1   = Z;
    Z       = -pinv(D)*(P - D)*Z + pinv(D);
    ztmp2   = Z;
    k       = k + 1; 
end
k
ztmp2;
inv(P);
ztmp2 - inv(P)

ptmp = zeros(n,1);
for i = 1:n
    for j = 1:n
        if i~=j
            ptmp(i)=ptmp(i)+P(i,j);
        end
    end
end
(ptmp<diag(P))'

U = triu(P,1);
L = tril(P,-1);
ztmp1 = inv(D);
ztmp2 = inv(D);
k = 1;
w = 1;
while (k==1)||( (norm( eig(ztmp1) - eig(ztmp2),Inf )>1e-10)&&(k<50) )
    ztmp1   = Z;
    Z       = w.*pinv(D+w.*L) - pinv(D+w.*L)*(w.*U + (w-1).*D)*Z;
    ztmp2   = Z;
    k       = k + 1; 
end
k
ztmp2;
inv(P);
ztmp2 - inv(P)