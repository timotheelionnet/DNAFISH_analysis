function y = half_pulse_shape( coeffs, xpts)
%ad hoc half-pulse shape consisting of a gaussian half followed by a plateau:

%first half for x<xmax
%y = (xpts < xmax).*( y0 + y1 * exp( -(xpts - xmax ).^2/(2*sig1^2)) ); 

%second half for x>=xmax
%y = y + (xpts >= xmax).*( y0 + y1 ); 

%coeffs are [ xmax, y0, y1, sig1]
%where 
    %xmax is the point of maximum (or minimum)
    %y0 is the -Inf baseline
    %y1 is the amplitude of the first half 
    %sig1 the sigma of the first half

%note that +Inf baseline (y0 + y1) 
%is different from -Inf baseline (y0)

xmax = coeffs(1);
y0 = coeffs(2);
y1 = coeffs(3);
sig1 = coeffs(4);


y = (xpts < xmax).*( y0 + y1 * exp( -(xpts - xmax ).^2/(2*sig1^2)) ); 

y = y + (xpts >= xmax).*( y0 + y1 ); 

y(isinf(y)) = realmax/1000;
end
