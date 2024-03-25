function I = multigaussian7(Coeffs,xpts)

A1 = Coeffs(1); A2 = Coeffs(2); A3 = Coeffs(3);  A4 = Coeffs(4); A5 = Coeffs(5); A6 = Coeffs(6); A7 = Coeffs(7);
I0 = Coeffs(8); s0 = Coeffs(9);

I = A1*exp ( -(xpts - 1*I0).^2/(2*1*s0^2));
I = I + A2*exp ( -(xpts - 2*I0).^2/(2*2*s0^2));
I = I + A3*exp ( -(xpts - 3*I0).^2/(2*3*s0^2));
I = I + A4*exp ( -(xpts - 4*I0).^2/(2*4*s0^2));
I = I + A5*exp ( -(xpts - 5*I0).^2/(2*5*s0^2));
I = I + A6*exp ( -(xpts - 6*I0).^2/(2*6*s0^2));
I = I + A7*exp ( -(xpts - 7*I0).^2/(2*7*s0^2));
end