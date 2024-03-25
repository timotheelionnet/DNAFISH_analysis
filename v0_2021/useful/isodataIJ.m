function threshold = isodataIJ( I )
    %public Object[] exec(ImagePlus imp, String myMethod, boolean noWhite, boolean noBlack, boolean doIwhite, boolean doIset, boolean doIlog , boolean doIstackHistogram ) {
    %I= uint8(I(:));
    threshold = -1;
    
    [data,xdata] = imhist(I);
    numel(data)
    maxbin = max(find(data>0));
    minbin = min(find(data>0));

    data2 = data(minbin:maxbin);    

    % Apply the selected algorithm
    if numel(data2) < 2
        threshold = 0;
    else 
        threshold = IJDefault(data2); % re-implemeted so we can ignore black/white and set the bright or dark objects
    end

    threshold = threshold + minbin; % add the offset of the histogram
    threshold = xdata(threshold);
end

function level = IJDefault(data ) 
    % Based on Original IJ implementation for compatibility.
    maxValue = numel(data);
    min = 1;
    while ((data(min)==0) && (min<maxValue))
        min= min+1;
    end
    
    max = maxValue;
    while ((data(max)==0) && (max>1))
        max = max -1;
    end
    min
    max
    if (min>=max) 
        level = numel(data)/2;
        return 
    end
    
    movingIndex = min;
    result = Inf;
    while ((movingIndex+1)<=result && movingIndex<max-1);
        
        sum1=0;
        sum2=0;
        sum3=0;
        sum4=0;
        for i=min:movingIndex
            sum1 = sum1 + i*data(i);
            sum2 = sum2 + data(i);
        end
        
        for i=movingIndex+1:max
            sum3 = sum3 + i*data(i);
            sum4 = sum4 + data(i);
        end
        result = (sum1/sum2 + sum3/sum4)/2.0;
        movingIndex = movingIndex+1;
    end 
    
    level = round(result);
    return 
end