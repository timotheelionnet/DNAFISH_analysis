function [Coeffsout, confintres,chi2] = fit_pulse_shape(xdata,ydata,sigmadata,guesscoeffs,lb,ub,exclude)
%fit data to ad hoc pulse shape consisting of 2 connected gaussian halves:
%first half for x<xmax
% y = (xpts < xmax).*( y0 + (ymax - y0) * exp( -(xpts - xmax ).^2/(2*sig1^2)) ); 
% 
%second half for x>=xmax
% y = y + (xpts >= xmax).*( yinf + (ymax - yinf) * exp( -(xpts - xmax ).^2/(2*sig2^2)) ); 

%coeffs are [ xmax, y0, y1, y2,sig1,sig2]
%where 
    %xmax is the point of maximum (or minimum)
    %y0 is the -Inf baseline
    %y1 is the amplitude of the first half
    %y2 is the amplitude of the second gaussian half
    %sig1 sig2 the sigmas of the respective halves

%note that +Inf baseline (y0 + y1 - y2) 
%is different from -Inf baseline (y1)

xdata = reshape(double(xdata),numel(xdata),1); 
ydata = reshape(double(ydata),numel(ydata),1); 
sigmadata = reshape(double(sigmadata),numel(sigmadata),1) ;

%guesscoeffs = [xmax_guess,y0_guess,ymax_guess,yinf_guess,sig1_guess,sig2_guess];

% options = optimset('TolX',.001,'MaxIter',1000);
% [Coeffsout, resnorm]=lsqcurvefit(@pulse_shape, guesscoeffs, xdata,ydata,lb,ub,options);

g = fittype( @(a,b,c,d,e,f, x) pulse_shape([a,b,c,d,e,f],x) );

sigmadata(sigmadata == 0) = min(sigmadata(sigmadata ~= 0))/100;


[fitres,gof] = fit(xdata,ydata,g,'Lower',lb,'Upper',ub,'StartPoint',guesscoeffs,'Weights',1./sigmadata.^2,'Robust','Bisquare',...
    'Display','off','Exclude',exclude);

Coeffsout = coeffvalues(fitres);
confintres = confint(fitres);
chi2 = sum((pulse_shape(Coeffsout,xdata) - ydata).^2./sigmadata.^2)/(numel(xdata)-numel(Coeffsout));
end

