function I = multigaussian2(Coeffs,xpts)

A1 = Coeffs(1); A2 = Coeffs(2); 
I10 = Coeffs(3); I20 = Coeffs(4); 
s10 = Coeffs(5); s20 = Coeffs(6); 

I = A1*exp ( -(xpts - I10).^2/(2*s10^2));
I = I + A2*exp ( -(xpts - I20).^2/(2*s20^2));


end