function [ mu,sigma] = bootstrap_distribution(data,itnum)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

m = zeros(itnum,1);
npts = numel(data);

for i=1:itnum
    Idx = ceil(rand(npts,1)*npts);
    m(i) = mean(data(Idx));
end

mu = mean(m);
sigma = std(m);

end

