function [ x0,x1,a,b,c,chi2,xfit,yfit ] = fit_square_bump_function(x,y,sigma,nbins)
%fits data: x, y, sigma (= error) all vectors with same size 
%to a square bump function:
    % y = a   if x<x0
    % y = b   if x>=x0 & x<x1
    % y = c   if x>=x1
    
    % possible x0 and x1 locations each take nbins possible discrete values
    % equally positioned between min(x) and max(x)

%OUTPUT: 
    %x0 = location of the first step
    %x1 = location of the second step
    %a  = value of the function x<x0
    %b  = value of the function x>=x0 & x<x1
    %c  = value of the function x>=x1
    %chi2 = chi squared value = sum_i{ ( y_i - f(x_i))^2 }
    %xfit = x
    %yfit = f(x)
    
%the function loops over x0 and x1 values, and computes the exact solution 
%to chi2 minimization for each x0; 
%eventually selects the best x0,x1 and corresponding a, b, c values.

x0 = (max(x) - min(x))*(0:nbins)/nbins + min(x);

chi2 = Inf;
amin = NaN;
bmin = NaN;
cmin = NaN;
x0min = NaN;
x1min = NaN;
xfit = x;
yfit = NaN;

if ~isempty(sigma(sigma ~= 0))
    sigma(sigma == 0) = min(sigma(sigma ~= 0))/100;
else
    sigma = ones(size(sigma));
end

for i=1:(numel(x0)-2)
    
    x1 = x0(i+1:end);
    
    for j=1:numel(x1)
        idx = x<x0(i);
        a = sum( y(idx)./sigma(idx).^2 ) / sum(1./sigma(idx).^2);
        tmpchi2 = sum( (y(idx) - a).^2./sigma(idx).^2 );

        idx = (x>=x0(i)) & (x<x1(j));
        b = sum( y(idx)./sigma(idx).^2 ) / sum(1./sigma(idx).^2);
        tmpchi2 = tmpchi2 + sum( (y(idx) - b).^2./sigma(idx).^2 );

        idx = (x>=x1(j));
        c = sum( y(idx)./sigma(idx).^2 ) / sum(1./sigma(idx).^2);
        tmpchi2 = tmpchi2 + sum( (y(idx) - c).^2./sigma(idx).^2 );
 
        if tmpchi2 < chi2 && c < b
            chi2 = tmpchi2;
            amin = a;
            bmin = b;
            cmin = c;
            x0min = x0(i);
            x1min = x1(j);
            xfit = x;
            yfit = (x<x0min)*a + ((x>=x0min)&(x<x1min))*b + (x>=x1min)*c;
        end
    end
end

a = amin;
b = bmin;
c = cmin;
x0 = x0min;
x1 = x1min;

end

