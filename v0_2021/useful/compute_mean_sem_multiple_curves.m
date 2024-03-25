function [ m,s ] = compute_mean_sem_multiple_curves(x,e)
%computes the mean and sem of many curves each from an independent
%experiment
% x is a npoints x ncurves array where each column is the data for one curve
% e is a npoints x ncurves array where each column is the s.e.m. for one curve  

% computes m, the npoints x 1 data curve as the weighed mean of all the data
% estimates s, the npoints x 1 s.e.m. as the weighed variance of the data
% around m divided by sqrt(n-1). (estimator seems to do a good job whether
% the error is mostly intra-sample or between samples).

if isempty(x) || isempty(e)
    m = [];
    s = [];
    disp('Mean/Sem of multiple datasets warning: found no data');
    return
    
elseif size(x,2) == 1
    m = x;
    s = e;
    disp('Mean/Sem of multiple datasets warning: found only one dataset');
    return
end

m = sum ( x./e.^2, 2 ) ./ sum( 1./e.^2, 2 );

dx = x - repmat(m,1,size(x,2)) ;

s = sqrt ( sum( dx.^2./e.^2 , 2 ) ./ sum( 1./e.^2, 2 )  );

s = s/sqrt(size(x,2)-1); 


end

