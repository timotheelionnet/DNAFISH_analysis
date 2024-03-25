function [I,I1,I2,I3,I4,I5] = get_multigaussian(xpts,Coeffs)
A1 = Coeffs(1); A2 = Coeffs(2); A3 = Coeffs(3);  A4 = Coeffs(4); A5 = Coeffs(5);  I0 = Coeffs(6); s0 = Coeffs(7);

I = multigaussian(Coeffs,xpts);
I1 = A1*exp ( -(xpts - 1*I0).^2/(2*1*s0^2));
I2 = A2*exp ( -(xpts - 2*I0).^2/(2*2*s0^2));
I3 = A3*exp ( -(xpts - 3*I0).^2/(2*3*s0^2));
I4 = A4*exp ( -(xpts - 4*I0).^2/(2*4*s0^2));
I5 = A5*exp ( -(xpts - 5*I0).^2/(2*5*s0^2));


end