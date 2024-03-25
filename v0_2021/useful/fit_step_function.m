function [ x0,a,b,chi2,xfit,yfit ] = fit_step_function(x,y,sigma,nbins)
%fits data: x, y, sigma (= error) all vectors with same size 
%to a step function:
    % y = a   if x<x0
    % y = b   if x>=x0

    % possible x0 locations take nbins possible discrete values
    % equally positioned between min(x) and max(x)

%OUTPUT: 
    %x0 = location of the step
    %a  = value of the function x<x0
    %b  = value of the function x>=x0
    %chi2 = chi squared value = sum_i{ ( y_i - f(x_i))^2 }
    %xfit = x
    %yfit = f(x)
    
%the function loops over x0 values, and computes the exact solution to chi2
%minimization for each x0; eventually selects the best x0 and corresponding
%a, b values.

x0 = (max(x) - min(x))*(0:nbins)/nbins + min(x);

chi2 = Inf;
amin = NaN;
bmin = NaN;
x0min = NaN;
xfit = x;
yfit = NaN;

if ~isempty(sigma(sigma ~= 0))
    sigma(sigma == 0) = min(sigma(sigma ~= 0))/100;
else
    sigma = ones(size(sigma));
end

for i=1:numel(x0)
    
    idx = x<x0(i);
    a = sum( y(idx)./sigma(idx).^2 ) / sum(1./sigma(idx).^2);
    
    tmpchi2 = sum( (y(idx) - a).^2./sigma(idx).^2 );
    
    idx = x>=x0(i);
    b = sum( y(idx)./sigma(idx).^2 ) / sum(1./sigma(idx).^2);
    
    tmpchi2 = tmpchi2 + sum( (y(idx) - b).^2./sigma(idx).^2 );

    if tmpchi2 < chi2
        chi2 = tmpchi2;
        amin = a;
        bmin = b;
        x0min = x0(i);
        xfit = x;
        yfit = (x<x0min)*a + (x>=x0min)*b;
    end
end

a = amin;
b = bmin;

x0 = x0min;
        
end

