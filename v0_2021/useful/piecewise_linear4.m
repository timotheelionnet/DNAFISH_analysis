function y = piecewise_linear4(coeffs,x)
%contiuous 4-piece linear function:
    % y = ax+b   if x<x0
    % y = cx+d   if x>=x0 & x<x1
    % y = ex+f   if x>=x1 & x<x2
    % y = gx+h   if x>=x2
    
%Coeffs are [x0,x1-x0,x2-x1,a,b,c,e,g];

%d,f and h are calculated using the continuity constraints.

x0 = coeffs(1);
dx = coeffs(2);
dx2 = coeffs(3);
a = coeffs(4);
b = coeffs(5);
c = coeffs(6);
e = coeffs(7);
g = coeffs(8);

x1 = x0 + dx;
x2 = x1 + dx2;
d = (a-c)*x0 + b;
f = (c-e)*x1+d;
h = (e-g)*x2+f;

y = double(x<x0) .* ( a*x + b ) + ...
    double(x>=x0 & x<x1) .* ( c*x + d ) + ...
    double(x>=x1 & x<x2) .* ( e*x + f ) + ...
    double(x>=x2) .* ( g*x + h );
end