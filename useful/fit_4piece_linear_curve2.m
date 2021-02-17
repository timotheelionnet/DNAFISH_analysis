function [ x0,x1,x2,a,b,c,d,e,f,g,h,chi2,xfit,yfit,Coeffsout ] = fit_4piece_linear_curve2(x,y,sigma,varargin)
%fits data: x, y, sigma (= error) all vectors with same size 
%to a contiuous 3 piece linear function:
    % y = ax+b   if x<x0
    % y = cx+d   if x>=x0 & x<x1
    % y = ex+f   if x>=x1 & x<x2
    % y = gx+h   if x>=x2
    %optional input: x0,x1,x2: impose guess values for x0, x1 and x2 as the result
    %depends heavily on initial conditions.

%OUTPUT: 
    %x0 = location of the first break
    %x1 = location of the second break
    %x2 = location of the third break
    %a  = value of the function ax+b x<x0
    %b  = value of the function ax+b x<x0
    %c  = value of the function cx+d x>=x0 & x<x1
    %d  = value of the function cx+d x>=x0 & x<x1
    %e  = value of the function ex+f x>=x1 & x<x2
    %f  = value of the function ex+f x>=x1 & x<x2
    %g  = value of the function gx+h x>=x2
    %h  = value of the function gx+h x>=x2
    %chi2 = chi squared value = sum_i{ ( y_i - f(x_i))^2 }
    %xfit = x
    %yfit = f(x)
    %Coeffsout: the list of actually fitted coefficients, can be used in
    %the function y = piecewise_linear4(coeffs,x)
        %x0 
        %x1-x0 
        %x2-x1 
        %a  
        %b  
        %c  
        %e 
        %g
    
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
idx3 = ~idx2 & ~idx1 & (x< x2);
if sum(idx3) == 0
    idx3(end-1) = 1;
end

idx4 = ~idx3 & ~idx2 & ~idx1 & (x< x2);
if sum(idx4) == 0
    idx4(end) = 1;
end

%note that the second coeff is the delta x1-x0, not x1;
%d and f are constrained bu continuity, so the fitted coeffs are really:
%x0; x1-x0; a; b; c; e;
guesscoeffs = [x0,x1-x0,x2-x1,...
    (max(y(idx1)) - min(y(idx1)))/(x0 - min(x)),mean(y(idx1)),...
    (max(y(idx2)) - min(y(idx2)))/(max(x) - min(x)),...
    (max(y(idx3)) - min(y(idx3)))/(max(x) - min(x)),...
    (max(y(idx4)) - min(y(idx4)))/(max(x) - min(x))];

mindeltax = min(diff(x))/10;
maxdeltay = (max(y) - min(y))/mindeltax;

lb = [min(x)+mindeltax,mindeltax,mindeltax,...
    -maxdeltay,-maxdeltay*max(x),...
    -maxdeltay,...
    -maxdeltay,...
    -maxdeltay];
    
ub = [max(x)-mindeltax,max(x)-min(x)-mindeltax,max(x)-min(x)-mindeltax,...
    maxdeltay,maxdeltay*max(x),...
    maxdeltay,...
    maxdeltay,...
    maxdeltay];

g = fittype( @(x0,dx,dx2,a,b,c,e,g,x) piecewise_linear4([x0,dx,dx2,a,b,c,e,g],x) );

sigma(sigma == 0) = min(sigma(sigma ~= 0))/100;


[fitres,gof] = fit(x,y,g,'Lower',lb,'Upper',ub,...
    'StartPoint',guesscoeffs,...
    'Weights',1./sigma.^2,...
    'Display','off',...
    'Algorithm','Trust-Region',...
    'TolX',0.01,...
    'TolFun',0.01);

Coeffsout = coeffvalues(fitres);
if numel(x)>numel(Coeffsout)
    confintres = confint(fitres);
end


x0 = Coeffsout(1);
x1 = Coeffsout(1) + Coeffsout(2);
x2 = x1 + Coeffsout(3);
a = Coeffsout(4);
b = Coeffsout(5);
c = Coeffsout(6);
d = (a-c)*x0+b;
e = Coeffsout(7);
f = (c-e)*x1+d;
g = Coeffsout(8);
h = (e-g)*x2+f;

chi2 = sum((piecewise_linear4(Coeffsout,x) - y).^2./sigma.^2);
xfit = x;
yfit = piecewise_linear4(Coeffsout,x);
end

