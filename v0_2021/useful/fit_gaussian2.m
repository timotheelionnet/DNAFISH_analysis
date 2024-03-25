function [Coeffsout, gof2] = fit_gaussian2(xdata,ydata)
%best if using x,y from a histogram
%Outputs Coeffs = [I0, mu, sigma]
%where ydata = I0 * exp [ -(xdata - mu).^2 / (2*sigma^2)    ]

%converting to column vectors
if isrow(xdata)
    xdata = xdata';
end

if isrow(ydata)
    ydata = ydata';
end

%estimate sigma using measured pdf integral
d = diff(xdata);
d = [d(1);d];
sig = sum(ydata.*d)/(max(ydata)*sqrt(2*pi));

%estimate x0 cumputing the mean of the measured pdf
mu = sum(xdata.*ydata.*d) / sum(ydata.*d);

%set weights proportional to values
weights = (ydata + (ydata == 0).*min(ydata(ydata ~= 0))).^2;  %assigning weight of 1 to points with value zero

s = fitoptions('Method','NonlinearLeastSquares',...
               'Display','off',...
               'Lower',[0,mu-2*sig,0],...
               'Upper',[3*max(ydata),mu+2*sig,3*sig],...
               'Startpoint',[max(ydata),mu,sig],...
               'Weights',weights);
         
f = fittype('a*exp ( -(x - b).^2/(2*c^2))','options',s);
[c2,gof2] = fit(xdata,ydata,f);  
Coeffsout = coeffvalues(c2);

end