function [y, umax]=findmaxF(rho,gamma)

%This function finds the maximum of the expression:
%Q(u)*Q(sqrt(1+rho^2)*gamma - rho*u) 
%using interval halving
y=[];
umax = NaN;

if rho<0
    warning('rho must be positive')
end
if rho==0
    y = normcdf(-gamma);
    return
end
if rho>1
    rho = 1/rho;
    
end


umax = sqrt(1+rho^2)*gamma/(1+rho);

umin = min(0,sqrt(1+rho^2)*gamma/rho - sqrt(2)/(sqrt(pi)*rho^2));

der_max = normpdf(umax)/normcdf(-umax) -... 
        rho*normpdf(sqrt(1+rho^2)*gamma-rho*umax)/normcdf(-sqrt(1+rho^2)*gamma+rho*umax);
% der_max = normpdf(umax)/normcdfFAST2(-umax) -... 
%         rho*normpdf(sqrt(1+rho^2)*gamma-rho*umax)/normcdfFAST2(-sqrt(1+rho^2)*gamma+rho*umax);
    
der_min = normpdf(umin)/normcdf(-umin) -... 
        rho*normpdf(sqrt(1+rho^2)*gamma-rho*umin)/normcdf(-sqrt(1+rho^2)*gamma+rho*umin);
% der_min = normpdf(umin)/normcdfFAST2(-umin) -... 
%         rho*normpdf(sqrt(1+rho^2)*gamma-rho*umin)/normcdfFAST2(-sqrt(1+rho^2)*gamma+rho*umin);

    
crit = 1e-12;
    u = .5*(umax+umin);

while abs(umax-umin)*max(abs(der_max),abs(der_min))>crit
    
    u = .5*(umax+umin);
    der =  normpdf(u)/normcdf(-u) -... 
        rho*normpdf(sqrt(1+rho^2)*gamma-rho*u)/normcdf(-sqrt(1+rho^2)*gamma+rho*u);
%     der =  normpdf(u)/normcdfFAST2(-u) -... 
%         rho*normpdf(sqrt(1+rho^2)*gamma-rho*u)/normcdfFAST2(-sqrt(1+rho^2)*gamma+rho*u);
    if der>0
        umax = u;
        der_max = der;
    else
        umin =u;
        der_min = u;
    end
end

maxb = u;
y = min(1,normcdf(-u)*normcdf(-sqrt(1+rho^2)*gamma+rho*u));
    
  

%plot function to be maximized
%   umax = sqrt(1+rho^2)*gamma/(1+rho);
% umin = min(0,sqrt(1+rho^2)*gamma/rho - sqrt(2)/(sqrt(pi)*rho^2));
%  u= umin:.1:umax; 
%   
%   y_f = normcdf(-u).*normcdf(-sqrt(1+rho^2)*gamma+rho*u);
% 
% figure
% semilogy(u,y_f)
% xlabel('u')
% ylabel('y')  
    
    
    
    
    
    