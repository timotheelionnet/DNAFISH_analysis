function [ x0,x1,a,b,c,d,chi2,xfit,yfit,Coeffsout ] = fit_3piece_linear_curve_flat_slope_flat2(x,y,sigma,varargin)
%fits data: x, y, sigma (= error) all vectors with same size 
%to a contiuous 3 piece linear function:
    % y = a   if x<x0
    % y = bx+c   if x>=x0 & x<x1
    % y = d   if x>=x1
    %optional input: x0,x1: impose guess values for x0 and x1 as the result
    %depends heavily on initial conditions.

%OUTPUT: 
    %x0 = location of the first step
    %x1 = location of the second step
    %a  = value of the function a x<x0
    %b  = value of the function bx+c x>=x0 & x<x1
    %c  = value of the function bx+c x>=x0 & x<x1
    %d  = value of the function d x>=x1
    %chi2 = chi squared value = sum_i{ ( y_i - f(x_i))^2 }
    %xfit = x
    %yfit = f(x)
    %Coeffsout: the list of actually fitted coefficients, can be used in
    %the function y = piecewise_linear3(coeffs,x):
        %[x0, x1-x0, 0, a, b, 0]   
    
[x,idx] = sort(x);
y = y(idx);
sigma = sigma(idx);

x = reshape(double(x),numel(x),1) ;
y = reshape(double(y),numel(y),1) ;
sigma = reshape(double(sigma),numel(sigma),1) ;

%separate dataset into 3 guess regions along x
if numel(varargin)>=2
    x0 = varargin{1};
    x1 = varargin{2};
    if x0<=min(x) || x0>=max(x)
        x0 = min(x) + (max(x) - min(x))/3;
    end
    
    if x0 >= x1 || x1 >=max(x)
        x1 = x0 + (max(x)-x0)/2;
    end
else  
    mindeltax = min(diff(x))/100;
    x0 = min(x) + (max(x)-min(x))/2 - mindeltax;
    x1 = x0 + mindeltax;
end

idx1 = x< x0;
if sum(idx1) == 0
    idx1(1) = 1;
end
idx2 = ~idx1 & (x< x1);
if sum(idx2) == 0
    if max(find(idx1))+1 <= numel(idx2)
        idx2(max(find(idx1))+1) = 1;
    else
        idx2(max(find(idx1))) = 1;
    end
end
idx3 = ~idx1 & ~idx2;
if sum(idx3) == 0
    idx3(end) = 1;
end
% x0 = x0
% x1 = x1
%note that the second coeff is the delta x1-x0, not x1;
%d and f are constrained bu continuity, so the fitted coeffs are really:
%x0; x1-x0; a; b; c; e;

guesscoeffs = [x0,x1-x0,...
    mean(y(idx1)),...
    (max(y) - min(y))/(x1-x0)];

mindeltax = min(diff(x))/100;
maxdeltay = (max(y) - min(y))/mindeltax;

lb = [min(x)+mindeltax,mindeltax,...
    -maxdeltay*max(x),...
    -maxdeltay];
    
ub = [max(x)-mindeltax,max(x)-mindeltax,...
    maxdeltay*max(x),...
    maxdeltay];

g = fittype( @(x0,dx,a,b,x) piecewise_linear3([x0,dx,0,a,b,0],x) );

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


x0 = Coeffsout(1);
x1 = Coeffsout(1) + Coeffsout(2);
a = Coeffsout(3);
b = Coeffsout(4);
c = piecewise_linear3([x0,x1-x0,0,a,b,0],x1) - b*x1;
d = piecewise_linear3([x0,x1-x0,0,a,b,0],x1);
xfit = x;
yfit = piecewise_linear3([x0,x1-x0,0,a,b,0],x);
chi2 = sum((piecewise_linear3([x0,x1-x0,0,a,b,0],x) - y).^2../sigma.^2);
Coeffsout = [x0,x1-x0,0,a,b,0];
end


