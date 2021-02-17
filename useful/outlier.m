function [dataclean,outliers,index] = outlier( y, alpha,k)
% Detection and Removal of Outliers in Data Sets
%    ( Rosner's many-outlier test)
%
%       index = outlier( y [, crit] )
%
% where  index = indices of outliers in the data
%        y     = data set (should be stationary)
%        crit  = detection criterion (default 2)
%
% Originally written by Bob Newell, February 1996
% Modified by Jaco de Groot, May 2006
% Bob Newell used a fixed value for lambda. This script calculates the
% critical values for lambda based on the equations in
% "Quality control of semi-continuous mobility size-fractionated particle number concentration data",
% Atmospheric Environment 38 (2004) 3341–3348, Rong Chun Yu,*, Hee Wen Teh, Peter A. Jaques, Constantinos Sioutas,
% John R. Froines)
%-----------------------------------------------------
%

y = y(:); 
n = length( y ); 
if nargin < 2 
    alpha = 0.05; 
else 
    if nargin < 3, 
        k = 1; 
    end 
end 
R = zeros( k+1, 1 ); 

%% sort deviations from the median 
ybar = median( y ); 
[ ys, is ] = sort( abs( y - ybar )); 

%% calculate statistics for up to k outliers 
for i = 0:k 
  yy = ys(1:n-i); 
  R(i+1) = abs( yy(n-i) - mean(yy) ) / std(yy); 
end

%% statistical test to find outliers 
index1 = []; index = []; 
lambda = zeros(k,1);
for i = 1:k 
    pcrit = 1 - alpha/(2*(n-i+1)); 
    t = tinv(pcrit, n-i-1); 
    
    lambda(i)=(n-i)*t./sqrt(((n-i-1+t^2)*(n-i+1))); 
    if R(i) > lambda 
        index=is(n-i+1:end); 
        index1 = [ index1 is(n-i+1) ]; 
    end 
end

idx = ones(1,numel(y));
idx(index1) = 0;
outliers = y(index1);
dataclean = y(logical(idx));

end


