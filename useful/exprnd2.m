function r = exprnd2(mu,varargin)

%EXPRND Random arrays from exponential distribution.
%   R = EXPRND(MU) returns an array of random numbers chosen from the
%   exponential distribution with mean parameter MU.  The size of R is
%   the size of MU.
%
%   R = EXPRND(MU,M,N,...) or R = EXPRND(MU,[M,N,...]) returns an
%   M-by-N-by-... array.
%
%   See also EXPCDF, EXPFIT, EXPINV, EXPLIKE, EXPPDF, EXPSTAT, RANDOM.

%   EXPRND uses the inversion method.

%   References:
%      [1]  Devroye, L. (1986) Non-Uniform Random Variate Generation, 
%           Springer-Verlag.

%   Copyright 1993-2009 The MathWorks, Inc. 
%   $Revision: 2.14.4.4 $  $Date: 2009/05/07 18:31:00 $

if nargin < 1
    error('stats:exprnd:TooFewInputs','Requires at least one input argument.');
end

sizeOut = size(zeros(varargin{1:end}));
    
% Return NaN for elements corresponding to illegal parameter values.
mu(mu < 0) = NaN;

% Generate uniform random values, and apply the exponential inverse CDF.
r = -mu .* log(rand(sizeOut)); % == expinv(u, mu)

clear('mu','err','sizeOut','tmp','argnum','numElements');

