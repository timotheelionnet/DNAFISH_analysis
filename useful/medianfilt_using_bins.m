function [ymed,xcenters,ystd,n,ysem,xmed] = medianfilt_using_bins(x,y,edges)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%reshape the edges as a line
edges = reshape(edges,1,numel(edges));

%centers of each bin
xcenters = (edges(1:end-1) + edges(2:end))/2;

ymed = zeros(numel(xcenters),1);
xmed = zeros(numel(xcenters),1);
ystd = zeros(numel(xcenters),1);
ysem = zeros(numel(xcenters),1);
n = zeros(numel(xcenters),1);

for i=1:numel(xcenters)
    idx = (x>edges(i)) & (x<= edges(i+1));
    ymed(i) = nanmedian(y(idx));
    xmed(i) = nanmedian(x(idx));
    n(i) = sum(idx);
    ystd(i) =  nanstd(y(idx));
    ysem(i) =  nanstd(y(idx))/sqrt(n(i));
end

end

