function [ a,b,c,d,x0,yfit,xfit] = fit_2piece_linear_bump(x,y,sigma,nbins)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

x0 = (max(x) - min(x))*(0:nbins)/nbins + min(x);

chi2 = Inf;
amin = NaN;
bmin = NaN;
cmin = NaN;
dmin = NaN;
x0min = NaN;
xfit = x;
yfit = NaN;


for i=1:numel(x0)
    
    tmpchi2 = 0;
    
    %find the coeffs of the first part of the curve
    idx = x<x0(i);
    sx = sum(x(idx)./sigma(idx).^2);
    sx2 = sum(x(idx).^2./sigma(idx).^2);
    sy = sum(y(idx)./sigma(idx).^2);
    sxy = sum(x(idx).*y(idx)./sigma(idx).^2);
    n = sum(1./sigma(idx).^2);    
    
    v = inv([sx2, sx; sx, n])*[sxy;sy];
    a = v(1);
    b = v(2);
    
    tmpchi2 = tmpchi2 + sum( ( y(idx) - a*x(idx) - b  ).^2./sigma(idx).^2  );
    
    %find the coeffs of the second part of the curve
    idx = x>=x0(i);
    sx = sum((x(idx) - x0(i))./sigma(idx).^2);
    sx2 = sum((x(idx) - x0(i)).^2./sigma(idx).^2);
    sxy = sum((x(idx) - x0(i)).*y(idx)./sigma(idx).^2);
    
    c = - ( sx*(a*x0(i) + b) - sxy )/sx2;
    d = (a-c)*x0(i) + b;
    
    tmpchi2 = tmpchi2 + sum( ( y(idx) - c*x(idx) - d  ).^2./sigma(idx).^2  );
    
    if tmpchi2 < chi2
        chi2 = tmpchi2;
        amin = a;
        bmin = b;
        cmin = c;
        dmin = d;
        x0min = x0(i);
        xfit = x;
        yfit = (x<x0min).*(a*x+b) + (x>=x0min).*(c*x+d);
    end
end

a = amin;
b = bmin;
c = cmin;
d = dmin;
x0 = x0min;
        

end

