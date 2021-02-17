function m = hist2d(xdata,ydata,xybins)

%INPUT
%xdata: x coordinates of spots (row or column vector)
%ydata: y coordinates of spots (row or column vector)
%xybins: coordinates of the centers of the bins along one axis
    %these coordinates will be used for both x and y axes.
    %the function assumes all bins are equally spaced

%gmask (optional, enter [] to ignore): normalize distances to the number of expected occurences based on a square image of size imsize    
%this option requires distances entered in pixel units
    
%imsize (optional, enter [] to ignore): normalize distances to the number of expected occurences based on a square image of size imsize    
%this option requires distances entered in pixel units

%spots outside of the interval [min(xybins),max(xybins)] are not plotted.


%OUTPUT
%m is the 2D histogram matrix (x vlues as rows, y values as columns)


%%
edges = hist_centers_to_edges(xybins,0);

%generate bin
binsize = xybins(2)-xybins(1);

%exclude data values outside of the grid
xymin = min(xybins ) - binsize/2;
xymax = max(xybins) + binsize/2;

idx = logical( (xdata>=xymin ).*(xdata<xymax ).*(ydata>=xymin ).*(ydata<xymax ));
xdata = xdata( idx );
ydata = ydata( idx );
clear('idx)');

xdata = reshape(xdata,numel(xdata),1);
ydata = reshape(ydata,numel(ydata),1);


%map points to bins
[~,~,subs] = histcounts([xdata,ydata], edges); 

m = accumarray(subs, 1, [numel(edges)-1, numel(edges)-1]);

end

