function y = piecewise_linear3(coeffs,x)
%contiuous 3-piece linear function:
    % y = ax+b   if x<x0
    % y = cx+d   if x>=x0 & x<x1
    % y = ex+f   if x>=x1
%Coeffs are [x0,x1-x0,a,b,c,e];

%d and f are calculated using the continuity constraints.

x0 = coeffs(1);
dx = coeffs(2);
a = coeffs(3);
b = coeffs(4);
c = coeffs(5);
e = coeffs(6);

x1 = x0 + dx;
d = (a-c)*x0 + b;
f = (c-e)*x1+d;

y = double(x<x0) .* ( a*x + b ) + ...
    double(x>=x0 & x<x1) .* ( c*x + d ) + ...
    double(x>=x1) .* ( e*x + f );
end