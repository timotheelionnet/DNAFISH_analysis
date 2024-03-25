function [yavg,ystd,n,xavg,xstd,xbin] = bin_ydata_using_xdata( xdata,ydata,xbin_ctrs,exclude_nan)
%combines x values into bins according to xbin_ctrs
%for each x bin, outputs the average y value of all datapoints that fall
%into the bin

%bin centeres should be sorted
z = diff(xbin_ctrs)<=0;
if sum(z) ~=0
    disp('bin centers should be sorted in ascending order');
    return
end
    
%compute boundaries between bins
bin_lims = (xbin_ctrs(2:end) + xbin_ctrs(1:numel(xbin_ctrs)-1))/2;


%assign xdata to each bin
xdata = reshape(xdata,numel(xdata),1);
bin_lims = reshape(bin_lims,1,numel(bin_lims));

idx_matrix = repmat(xdata,1,numel(bin_lims)) <= repmat(bin_lims,numel(xdata),1);
idx_matrix = sum(idx_matrix,2);
idx_matrix = numel(bin_lims)-idx_matrix+1;

%compile values in bins
for i=1:numel(xbin_ctrs)
    x = idx_matrix == i;
    if exclude_nan
        yavg(i) = nanmean(ydata(x));
        ystd(i) = nanstd(ydata(x));
    else
        yavg(i) = mean(ydata(x));
        ystd(i) = std(ydata(x));
    end
    
    n(i) = sum(x);
    xavg(i) = mean(xdata(x));
    xstd(i) = std(xdata(x));
       
end
xbin = xbin_ctrs;
end

