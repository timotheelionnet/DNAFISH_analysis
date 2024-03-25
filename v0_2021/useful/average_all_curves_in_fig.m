function res = average_all_curves_in_fig(fh,xbin)

%% collect all curves

%panels
axesObjs = get(fh, 'Children');

%loop through panels
for j=1:numel(axesObjs)
    %get data X and Y
    dataObjs = get(axesObjs(j), 'Children');
    xdata = get(dataObjs, 'XData');
    ydata = get(dataObjs, 'YData');
    
    %extract each curve separately
    res.xmax(j) = -Inf;
    res.xmin(j) = Inf;
    
    for k= 1:size(xdata,1)
        if iscell(xdata)
            res.xdata{j,k} = xdata{k};
            res.ydata{j,k} = ydata{k};
        else
            res.xdata{j,k} = xdata(k,:);
            res.ydata{j,k} = ydata(k,:);
        end
        
        if min(res.xdata{j,k}) < res.xmin(j)
            res.xmin(j) = min(res.xdata{j,k});
        end
        if max(res.xdata{j,k}) > res.xmax(j)
            res.xmax(j) = max(res.xdata{j,k});
        end
    end
end

%compute interval for each panel
for j=1:numel(axesObjs)
    j
    xbmin = xbin*floor(res.xmin(j)/xbin)
    xbmax = xbin*ceil(res.xmax(j)/xbin)
    res.xbins{j} = xbmin:xbin:xbmax;
    
    res.xavg{j} = NaN*zeros(1,numel(res.xbins{j})-1);
    res.yavg{j} = NaN*zeros(1,numel(res.xbins{j})-1);
    res.xstd{j} = NaN*zeros(1,numel(res.xbins{j})-1);
    res.ystd{j} = NaN*zeros(1,numel(res.xbins{j})-1);
    res.xsem{j} = NaN*zeros(1,numel(res.xbins{j})-1);
    res.ysem{j} = NaN*zeros(1,numel(res.xbins{j})-1);
    res.ncurves{j} = zeros(1,numel(res.xbins{j})-1);
    
    for ib= 1 : (numel(res.xbins{j})-1)
        
        for k=1:size(res.xdata,2)
            
            if ~isempty(res.xdata{j,k})
                xIdx = logical( (res.xdata{j,k}>=res.xbins{j}(ib)) .* (res.xdata{j,k}<res.xbins{j}(ib+1)) );
                res.ncurves{j}(ib) = res.ncurves{j}(ib) + sum(double(xIdx));
                if isnan(res.xavg{j}(ib))
                    res.xavg{j}(ib) = sum(res.xdata{j,k}(xIdx));
                    res.yavg{j}(ib) = sum(res.ydata{j,k}(xIdx));
                    res.xstd{j}(ib) = sum(res.xdata{j,k}(xIdx).^2);
                    res.ystd{j}(ib) = sum(res.ydata{j,k}(xIdx).^2);
                else
                    res.xavg{j}(ib) = res.xavg{j}(ib)+ sum(res.xdata{j,k}(xIdx));
                    res.yavg{j}(ib) = res.yavg{j}(ib)+ sum(res.ydata{j,k}(xIdx));
                    res.xstd{j}(ib) = res.xstd{j}(ib)+ sum(res.xdata{j,k}(xIdx).^2);
                    res.ystd{j}(ib) = res.ystd{j}(ib)+ sum(res.ydata{j,k}(xIdx).^2);
                end
            end
        end
        
    end
    res.xavg{j} = res.xavg{j}./res.ncurves{j};
    res.yavg{j} = res.yavg{j}./res.ncurves{j};
   
    res.xstd{j} = sqrt(res.xstd{j}./res.ncurves{j} - res.xavg{j}.^2);
    res.ystd{j} = sqrt(res.ystd{j}./res.ncurves{j} - res.yavg{j}.^2);
   
    res.xsem{j} = res.xstd{j}./sqrt(res.ncurves{j});
    res.ysem{j} = res.ystd{j}./sqrt(res.ncurves{j});
   
    axes(axesObjs(j));
    hold(axesObjs(j),'on');
    errorbar(res.xavg{j},res.yavg{j},res.ysem{j});
end






end
