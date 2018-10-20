function dF = dfdx(f,k,x,u,w,v,n);
% dF = dfdx(f,k,x,u,w,v,n);
% f   Function to be differentiated 
% k   time instant
% u   Input
% w   Process Noise
% v   Measurement Noise
% n   Order of the System

% Finite diffence step
hx = 1.0e-8;

% Gradient of objective function
fx = f(k,x,u,w,v);
fx1plush =cell(n,1);
dF = [];
for i = 1:n
    tmp         = zeros(n,1);
    tmp(i)      = hx;
    fx1plush{i} = f(k,x+hx,u,w,v);
    dfdx1       = (fx1plush{i} - fx)/hx;
    dF          = [dF,dfdx1];
end 
% fx1plush = f(k,[x(1)+hx x(2)]',u,w,v);
% fx2plush = f(k,[x(1), x(2)+hx]',u,w,v);
% dfdx1 = (fx1plush - fx)/hx;
% dfdx2 = (fx2plush - fx)/hx;
% dF = [dfdx1 dfdx2];
end 