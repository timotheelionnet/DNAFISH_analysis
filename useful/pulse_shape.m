function y = pulse_shape( coeffs, xpts)
%ad hoc pulse shape consisting of 2 connected gaussian halves:
%first half for x<xmax
%y = (xpts < xmax).*( y0 + y1 * exp( -(xpts - xmax ).^2/(2*sig1^2)) ); 

%second half for x>=xmax
%y = y + (xpts >= xmax).*( (y0 + y1 - y2) + y2 * exp( -(xpts - xmax ).^2/(2*sig2^2)) ); 

%coeffs are [ xmax, y0, y1, y2,sig1,sig2]
%where 
    %xmax is the point of maximum (or minimum)
    %y0 is the -Inf baseline
    %y1 is the amplitude of the first half
    %y2 is the amplitude of the second gaussian half
    %sig1 sig2 the sigmas of the respective halves

%note that +Inf baseline (y0 + y1 - y2) 
%is different from -Inf baseline (y0)

xmax = coeffs(1);
y0 = coeffs(2);
y1 = coeffs(3);
y2 = coeffs(4);
sig1 = coeffs(5);
sig2 = coeffs(6);

y = (xpts < xmax).*( y0 + y1 * exp( -(xpts - xmax ).^2/(2*sig1^2)) ); 

y = y + (xpts >= xmax).*( (y0 + y1 - y2) + y2 * exp( -(xpts - xmax ).^2/(2*sig2^2)) ); 

y(isinf(y)) = realmax/1000;
end

