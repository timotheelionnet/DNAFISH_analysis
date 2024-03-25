function [ x0,x1,x2,a,b,c,d,e,f,chi2,xfit,yfit,Coeffsout ] = fit_4piece_linear_curve_flat_slope_slope_flat2(x,y,sigma,varargin)
%fits data: x, y, sigma (= error) all vectors with same size 
%to a contiuous 3 piece linear function:
    % y = a   if x<x0
    % y = bx+c   if x>=x0 & x<x1
    % y = dx+e   if x>=x1 & x<x2
    % y = f   if x>=x2
    %optional input: x0,x1,x2: impose guess values for x0, x1 and x2 as the result
    %depends heavily on initial conditions.

%OUTPUT: 
    %x0 = location of the first break
    %x1 = location of the second break
    %x2 = location of the third break
    %a  = value of the function a x<x0
    %b  = value of the function bx+c x>=x0 & x<x1
    %c  = value of the function bx+c x>=x0 & x<x1
    %d  = value of the function dx+e x>=x1 & x<x2
    %e  = value of the function dx+e x>=x1 & x<x2
    %f  = value of the function f x>=x2
    %f  = value of the function f x>=x2
    %chi2 = chi squared value = sum_i{ ( y_i - f(x_i))^2 }
    %xfit = x
    %yfit = f(x)
    %Coeffsout: the list of actually fitted coefficients, can be used in
    %the function y = piecewise_linear4(coeffs,x)
        %[x0 , x1-x0 , x2-x1 , 0 , a , c , e , 0];
    
x = reshape(double(x),numel(x),1) ;
y = reshape(double(y),numel(y),1) ;
sigma = reshape(double(sigma),numel(sigma),1) ;

%separate dataset into 4 guess regions along x
if numel(varargin)>=3
    x0 = varargin{1};
    x1 = varargin{2};
    x2 = varargin{3};
    
    if x0<=min(x) || x0>=max(x)
        x0 = min(x) + (max(x) - min(x))/4;
    end
    
    if x0 >= x1 || x1 >=max(x)
        x1 = x0 + (max(x)-x0)/3;
    end
    
    if x1 >= x2 || x2 >=max(x)
        x2 = x1 + (max(x)-x1)/3;
    end
else  
    mindeltax = min(diff(x))/100;
    x0 = min(x) + (max(x)-min(x))/4;
    x1 = min(x) + 2*(max(x)-min(x))/4;
    x2 = min(x) + 3*(max(x)-min(x))/4;
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

%note that the second coeff is the delta x1-x0, not x1;
%d and f are constrained bu continuity, so the fitted coeffs are really:
%x0; x1-x0; a; b; c; e;
guesscoeffs = [x0,x1-x0,x2-x1,...
    mean(y(idx1)),...
    (max(y(idx2)) - min(y(idx2)))/(x1-x0),...
    (max(y(idx3)) - min(y(idx3)))/(x2-x1)];

mindeltax = min(abs(diff(x)))/10;
maxdeltay = (max(y) - min(y))/mindeltax;

lb = [min(x)+mindeltax,mindeltax,mindeltax,...
    -maxdeltay*max(x),...
    -maxdeltay,...
    -maxdeltay];
    
ub = [max(x)-mindeltax,max(x)-min(x)-mindeltax,max(x)-min(x)-mindeltax,...
    maxdeltay*max(x),...
    maxdeltay,...
    maxdeltay];

g = fittype( @(x0,dx,dx2,a,b,d,x) piecewise_linear4([x0,dx,dx2,0,a,b,d,0],x) );

sigma(sigma == 0) = min(sigma(sigma ~= 0))/100;


[fitres,gof] = fit(x,y,g,'Lower',lb,'Upper',ub,...
    'StartPoint',guesscoeffs,...
    'Weights',1./sigma.^2,'Robust','Bisquare',...
    'Display','off');

Coeffsout = coeffvalues(fitres);
if numel(x)>numel(Coeffsout)
    confintres = confint(fitres);
end


x0 = Coeffsout(1);
x1 = Coeffsout(1) + Coeffsout(2);
x2 = x1 + Coeffsout(3);
a = Coeffsout(4);
b = Coeffsout(5);
d = Coeffsout(6);
c = piecewise_linear4([x0,x1-x0,x2-x1,0,a,b,d,0],x0) - b*x0;
e = piecewise_linear4([x0,x1-x0,x2-x1,0,a,b,d,0],x1) - d*x1;
f = piecewise_linear4([x0,x1-x0,x2-x1,0,a,b,d,0],x2) ;

%these coeffs can be fed into piecewise_linear4
Coeffsout = [x0,x1-x0,x2-x1,0,a,b,d,0];
xfit = x;
yfit = piecewise_linear4(Coeffsout,x);

chi2 = sum((piecewise_linear4(Coeffsout,x) - y).^2./sigma.^2);
end

