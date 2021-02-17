function [edges,xc] = hist_centers_to_edges(xc,enforce_positive_values)
%Converts positions of bin centers xc into edges values
%e.g. [0.5,1.5,2.5]
%yields [0,1,2,3]

%edges are set at the midpoint between each pair of adjacent bin centers

%edges of first and last bin are set using the same bin width as the
%neighboring bin

%enforce_positive_values returns only positive bins


%edages are the edges of the bins
%xc is the potentially modified values of the corresponding bin centers

%format as line vector
xc =reshape(xc,1,numel(xc));

%remove negatives
if enforce_positive_values
    xc(xc<0) = 0;
    xc = unique(xc);
end

%sort
xc = sort(xc,'ascend');

bin_boundaries = (xc(1:end-1) + xc(2:end))/2;
bmin = [2*xc(1) - bin_boundaries(1),bin_boundaries];
if enforce_positive_values
    bmin = max(0,bmin);
end

bmax = [bin_boundaries, 2*xc(end) - bin_boundaries(end) ];
edges = [bmin, bmax(end)];

end

