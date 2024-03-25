function [ x0,x1,a,b,c,d,e,f,chi2,xfit,yfit,Coeffsout ] = fit_3piece_linear_curve(x,y,sigma)
%fits data: x, y, sigma (= error) all vectors with same size 
%to a contiuous 3 piece linear function:
    % y = ax+b   if x<x0
    % y = cx+d   if x>=x0 & x<x1
    % y = ex+f   if x>=x1

%OUTPUT: 
    %x0 = location of the first step
    %x1 = location of the second step
    %a  = value of the function ax+b x<x0
    %b  = value of the function ax+b x<x0
    %c  = value of the function cx+d x>=x0 & x<x1
    %d  = value of the function cx+d x>=x0 & x<x1
    %e  = value of the function ex+f x>=x1
    %f  = value of the function ex+f x>=x1
    %chi2 = chi squared value = sum_i{ ( y_i - f(x_i))^2 }
    %xfit = x
    %yfit = f(x)
    %Coeffsout: the list of actually fitted coefficients, can be used in
    %the function y = piecewise_linear3(coeffs,x)
        %x0 
        %x1-x0 
        %a  
        %b  
        %c  
        %e  
    
x = reshape(double(x),numel(x),1) ;
y = reshape(double(y),numel(y),1) ;
sigma = reshape(double(sigma),numel(sigma),1) ;

%separate dataset into 3 equal regions along x
x0 = min(x) + (max(x)-min(x))/3;
x1 = min(x) + 2*(max(x)-min(x))/3;

idx1 = x< x0;
idx2 = ~idx1 & (x< x1);
idx3 = ~idx1 & ~idx2;

%note that the second coeff is the delta x1-x0, not x1;
%d and f are constrained bu continuity, so the fitted coeffs are really:
%x0; x1-x0; a; b; c; e;
guesscoeffs = [x0,x1-x0,...
    (max(y(idx1)) - min(y(idx1)))/(x0 - min(x)),mean(y(idx1)),...
    (max(y(idx2)) - min(y(idx2)))/(x1-x0),...
    (max(y(idx3)) - min(y(idx3)))/(max(x)-x1)];

mindeltax = min(diff(x))/10;
maxdeltay = (max(y) - min(y))/mindeltax;

lb = [min(x)+mindeltax,mindeltax,...
    -maxdeltay,-maxdeltay*max(x),...
    -maxdeltay,...
    -maxdeltay];
    
ub = [max(x)-mindeltax,max(x)-min(x)-mindeltax,...
    maxdeltay,maxdeltay*max(x),...
    maxdeltay,...
    maxdeltay];

g = fittype( @(x0,dx,a,b,c,e,x) piecewise_linear3([x0,dx,a,b,c,e],x) );

sigma(sigma == 0) = min(sigma(sigma ~= 0))/100;


[fitres,gof] = fit(x,y,g,'Lower',lb,'Upper',ub,...
    'StartPoint',guesscoeffs,...
    'Weights',1./sigma.^2,...
    'Display','off',...
    'Algorithm','Trust-Region',...
    'TolX',0.01,...
    'TolFun',0.01);

Coeffsout = coeffvalues(fitres);
confintres = confint(fitres);
chi2 = sum((piecewise_linear3(Coeffsout,x) - y).^2../sigma.^2);

x0 = Coeffsout(1);
x1 = Coeffsout(1) + Coeffsout(2);
a = Coeffsout(3);
b = Coeffsout(4);
c = Coeffsout(5);
d = (a-c)*x0+b;
e = Coeffsout(6);
f = (c-e)*x1+d;
xfit = x;
yfit = piecewise_linear3(Coeffsout,x);
end

