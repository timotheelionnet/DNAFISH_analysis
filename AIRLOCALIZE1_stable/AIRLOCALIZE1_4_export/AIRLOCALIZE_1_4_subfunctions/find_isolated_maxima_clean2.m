function maxima = find_isolated_maxima_clean2(stack,threshInt,ROIsize)

    ROIsize = round(ROIsize);
    if ROIsize >0
        %find local maxima within ROI pixels in each dimension
        if ndims(stack) == 3
            %generate local window for dilation operation
            dilwin = ones(2*ROIsize+1,2*ROIsize+1,2*ROIsize+1);
            dilwin(ROIsize+1,ROIsize+1,ROIsize+1) = 0;
        elseif ndims(stack) == 2
            %generate local window for dilation operation
            dilwin = ones(2*ROIsize+1,2*ROIsize+1);
            dilwin(ROIsize+1,ROIsize+1) = 0;
        end
        localmax = stack>= imdilate(stack,dilwin);
    else
        localmax = ones(size(stack));
    end
    %enforce that local maxima intensity be above threshold 
    localmax = localmax.*(stack > threshInt);
    %store maxima as a list of coordinates / intensity
    x = find(localmax);
    clear('localmax');
    Int = stack(x);
    if ndims(stack) == 3
        [x,y,z] = ind2sub(size(stack),x);
        maxima = [x,y,z,Int];
    else
        [x,y] = ind2sub(size(stack),x);
        maxima = [x,y,Int];
    end
    
       
%%%%old function (inefficient code)
%    pts = find_points_above_threshold2(stack,threshInt);
%    maxima = remove_non_localmax_pixels2(stack,pts,ROIsize);
%%%%%    
    
%     if (ROIsize ~=0)
%         removes a point if there is another one with a higher intensity
%         within 'ROIsize' pixels in any direction (faster)
%         maxima = remove_non_localmax_pixels2(stack,pts,ROIsize);
%     else
%         for each cluster of neighbouring pixels, selects the one with the
%         highest intensity (slower) ABANDONNED
%         maxima = remove_connected_pixels2(stack,pts,threshInt);
%     end
    
end