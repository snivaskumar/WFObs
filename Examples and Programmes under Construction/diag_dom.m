function DP = diag_dom(A);
% DP = diag_dom(A)
% DP    Percentage of dominance

n = length(A);
D = diag(A);
R = A - diag(D);
RR = abs(R);
S = sum(RR');
S = S';
DD = double(D>=S);  
DS = sum(DD);
DP = (DS/n)*100;
end

% clear all
% close all
% clc
% 
% n = 6;
% L = 1;
% % A = rand(n,n);
% A = rand(n,n) + n*eye(n,n);
% A = (A*A')/norm(A*A');
% % bandwidth(A)
% 
% 
% 
% hr = 4;
% x       = cell(hr,1);
% % x{1}    = [1,2,3]';
% % x{2}    = [1,4]';
% % x{3}    = [4,5]';
% % x{4}    = [2,3,6]';
% 
% x{1}    = [1,2,3]';
% x{2}    = [3,4]';
% x{3}    = [4,5]';
% x{4}    = [5,6]';
% 
% % hr = 3;
% % x       = cell(hr,1);
% % x{1}    = [1,2,3]';
% % x{2}    = [2,3,4]';
% % x{3}    = [4,5]';
% 
% Z = inv(A);
% Z1 = inv(A(x{1},x{1}));
% Z2 = inv(A(x{2},x{2}));
% Z3 = inv(A(x{3},x{3}));
% Z(x{1},x{1})
% Z1
% Z(x{2},x{2})
% Z2
% Z(x{3},x{3})
% Z3
% 
% Z4 = inv(A(x{4},x{4}))
% Z(x{4},x{4})

% P = A;
% P
% P1 = P(x{1},x{1})