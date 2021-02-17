function [Coeffsout, resnorm] = fit_multigaussian(xdata,ydata,A1guess,A2guess,A3guess,A4guess,A5guess,I0guess,s0guess)
xdata = double(xdata); ydata = double(ydata);

Coeffs = [A1guess,A2guess,A3guess,A4guess,A5guess,I0guess,s0guess];
lb = [0,0,0,0,0,0.1,10];
ub = [1e6,1e6,1e6,1e6,1e6,1000,10000];

options = optimset('TolX',.001,'MaxIter',1000);
[Coeffsout, resnorm]=lsqcurvefit(@multigaussian, Coeffs, xdata,ydata,lb,ub,options);

end