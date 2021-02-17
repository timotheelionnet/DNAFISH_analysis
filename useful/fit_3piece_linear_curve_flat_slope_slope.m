function [ x0,x1,a,b,c,d,e,chi2,xfit,yfit,Coeffsout ] = fit_3piece_linear_curve_flat_slope_slope(x,y,sigma)
%fits data: x, y, sigma (= error) all vectors with same size 
%to a contiuous 3 piece linear function:
    % y = a   if x<x0
    % y = bx+c   if x>=x0 & x<x1
    % y = dx+e   if x>=x1

%OUTPUT: 
    %x0 = location of the first break
    %x1 = location of the second break
    %a  = value of the function a x<x0
    %b  = value of the function bx+c x>=x0 & x<x1
    %c  = value of the function bx+c x>=x0 & x<x1
    %d  = value of the function dx+e x>=x1
    %e  = value of the function dx+e x>=x1
    %chi2 = chi squared value = sum_i{ ( y_i - f(x_i))^2 }
    %xfit = x
    %yfit = f(x)
    %Coeffsout: the list of actually fitted coefficients, can be used in
    %the function y = piecewise_linear3(coeffs,x)
        %[x0, x1-x0, 0, a, b, d]
    
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
    mean(y(idx1)),...
    (max(y(idx2)) - min(y(idx2)))/(x1-x0),...
    (max(y(idx2)) - min(y(idx2)))/(x1-x0)];

mindeltax = min(diff(x))/10;
maxdeltay = (max(y) - min(y))/mindeltax;

lb = [min(x)+mindeltax,mindeltax,...
    -maxdeltay,...
    -maxdeltay*max(x),...
    -maxdeltay*max(x)];
    
ub = [max(x)-mindeltax,max(x)-min(x)-mindeltax,...
    maxdeltay,...
    maxdeltay*max(x),...
    maxdeltay*max(x)];

g = fittype( @(x0,dx,a,b,d,x) piecewise_linear3([x0,dx,0,a,b,d],x) );

sigma(sigma == 0) = min(sigma(sigma ~= 0))/100;


[fitres,gof] = fit(x,y,g,'Lower',lb,'Upper',ub,...
    'StartPoint',guesscoeffs,...
    'Weights',1./sigma.^2,'Robust','Bisquare',...
    'Display','off');

Coeffsout = coeffvalues(fitres);
confintres = confint(fitres);


x0 = Coeffsout(1);
x1 = Coeffsout(1) + Coeffsout(2);
a = Coeffsout(3);
b = Coeffsout(4);
d = Coeffsout(5);
c = piecewise_linear3([x0,x1-x0,0,a,b,d],x1) - b*x1;
e = piecewise_linear3([x0,x1-x0,0,a,b,d],x1) - d*x1;
xfit = x;
yfit = piecewise_linear3([x0,x1-x0,0,a,b,d],x);
chi2 = sum((piecewise_linear3([x0,x1-x0,0,a,b,d],x) - y).^2../sigma.^2);


Coeffsout = [x0,x1-x0,0,a,b,d];
end
