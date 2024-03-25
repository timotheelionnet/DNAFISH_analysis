function params = predict_threshold(params)
    %% retrieve data
    params = retrieve_smoothed_data(params);

    
    %% generate a set of threshold values spanning the image intensity
    %range
    npts = 100;
    Imax = double(max(params.smooth(:)));
    Imin = double(min(params.smooth(:)));
    thr = (Imax - Imin)*(0:npts-1)/(npts-1) + Imin;

    
    %% compute the number of objects above threshold for each value
    nobj = zeros(npts,1);
    for i =1:npts
        %disp(num2str(i));
        cc = bwconncomp(params.smooth > thr(i),18);
        nobj(i,1) = cc.NumObjects;  
        clear('cc');
    end
    
    res.thr = thr;
    
    %% finding the minimum threshold
    
    %index of the threshold value that maximizes the number of objects
    %detected (corresponds to the median background)
    [~,imax] = max(nobj);
    
    %index of the minimum threshold value thats picks up more than one
    %object (low end of the background)
    imin = find(nobj>1,1);
    params.thresh.thr_min = thr(imin);
    
    %"optimal index" : assuming the number of background objects grows from
    %imin to imax and then decreases symetrically, the optimal index
    %ensures that no background spots are detected.
    iopt = 2*imax - imin;
    params.thresh.thr_opt = thr(iopt);
    
    
end


