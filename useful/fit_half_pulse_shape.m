function [Coeffsout, confintres,chi2] = fit_half_pulse_shape(xdata,ydata,sigmadata,guesscoeffs,lb,ub,exclude)
%fit to ad hoc half-pulse shape consisting of a gaussian half followed by a plateau:

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

xdata = reshape(double(xdata),numel(xdata),1) ;
ydata = reshape(double(ydata),numel(ydata),1) ;
sigmadata = reshape(double(sigmadata),numel(sigmadata),1) ;

%guesscoeffs = [xmax_guess,y0_guess,ymax_guess,yinf_guess,sig1_guess,sig2_guess];

% options = optimset('TolX',.001,'MaxIter',1000);
% [Coeffsout, resnorm]=lsqcurvefit(@pulse_shape, guesscoeffs, xdata,ydata,lb,ub,options);

g = fittype( @(a,b,c,d, x) half_pulse_shape([a,b,c,d],x) );

if ~isempty(sigmadata(sigmadata ~= 0))
    sigmadata(sigmadata == 0) = min(sigmadata(sigmadata ~= 0))/100;
else
    sigmadata = ones(size(sigmadata));
end

[fitres,gof] = fit(xdata,ydata,g,'Lower',lb,'Upper',ub,'StartPoint',guesscoeffs,'Weights',1./sigmadata.^2,'Robust','Bisquare',...
    'Display','off','Exclude',exclude);

Coeffsout = coeffvalues(fitres);
confintres = confint(fitres);
chi2 = sum((half_pulse_shape(Coeffsout,xdata) - ydata).^2./sigmadata.^2);
end
