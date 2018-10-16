clear all
clc

s = 2*[1 1];
x = randn(5000,1);
x1 = normrnd(s(1).*x,1);
x2 = normrnd(s(2).*x,1);
x = [x1 x2];

[r_ellipse,X0,Y0,x,mu,var] = error_ellipse(x,1);

%figure,plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,'k-','linewidth',2),hold on,
%plot(mu(1),mu(2),'m+','linewidth',2)
close all;
%% Actual
y = (sin(x)).^2;    % Actual Function
% y = sin(4*x) - (x.^3)/10 + (2*x)/5;
[yr_ellipse,yX0,yY0,ydata,ymu] = error_ellipse(y,1);

plot(yr_ellipse(:,1) + yX0,yr_ellipse(:,2) + yY0,'k-','linewidth',2), hold on,
plot(ymu(1),ymu(2),'k+','linewidth',2)

%% EKF
A = diag(2*((sin(mu)).^1).*cos(mu));
%A = diag(2*((sin(mu)).^1).*cos(mu)) + [0 2*((sin(mu(2))).^1).*cos(mu(1)); 2*((sin(mu(1))).^1).*cos(mu(2)) 0];
% A = diag(4*cos(4*mu) - 3*(mu.^2)/10 -2/5);
ey = A*x';          % 1st order Linear Function
% f = (sin(mu)).^2;
% Df = 2*((sin(mu)).^1).*cos(mu);
% ey = f' + Df*(x - mu)';
[er_ellipse,eX0,eY0,edata,emu] = error_ellipse(ey',1);
plot(er_ellipse(:,1) + eX0,er_ellipse(:,2) + eY0,'b-','linewidth',2),
plot(emu(1),emu(2),'b+','linewidth',2)

%% UKF

mu = mean(x);
var = cov(x);
[rx cx] = size(x);
n = cx;
%w = (1/5)*ones(1,5);
alpha   = .83;
beta    = 2;
k       = 0.01;
% alpha   = 1;
% beta    = 0;
% k       = .8;
lambda  = (alpha^2)*(n + k) - n;
%ux = [mu; repmat(mu,1,n) + sqrt((2/(1-w))*var); repmat(mu,1,n) - sqrt((2/(1-w))*var)];
for i = 1:(2*n+1)
    if i == 1
        Wm(i)    = lambda/(n + lambda);
        Wc(i)    = Wm + (1 - alpha^2 + beta);
    else 
        Wm(i)    = 0.5/(n + lambda);
        Wc(i)    = Wm(i);
    end
end
% ux              = repmat(mu',1,2*n+1); 
% Uscented_devs   = sqrt(n+lambda)*sqrt(var)
% ux(:,2:n+1)     = ux(:,2:n+1) + Uscented_devs;
% ux(:,n+2:end)   = ux(:,n+2:end) - Uscented_devs;
ux = [mu; mu + sqrt(n+lambda)*sqrt(var); mu - sqrt(n+lambda)*sqrt(var)];
uy              = (sin(ux)).^2;  % UT
% uy = sin(4*ux) - (ux.^3)/10 + (2*ux)/5;

[uyr_ellipse,uyX0,uyY0,uydata,uymu] = error_ellipse(uy,0,Wm,Wc);
plot(uyr_ellipse(:,1) + uyX0,uyr_ellipse(:,2) + uyY0,'r-','linewidth',2),
plot(uymu(1),uymu(2),'r+','linewidth',2)