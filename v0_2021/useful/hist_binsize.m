function [n,x] = hist_binsize(arr,binsize)
%quick and dirty aproximate of histogram where you decide of the bin width

x1 = min(arr);
x2 = max(arr);

nbins = ceil((x2-x1)/binsize);

[n,x] = hist(arr,nbins);


end