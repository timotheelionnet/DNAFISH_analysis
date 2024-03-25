function [stack2 npix_above_thresh] = change_overlay_channel_clean(smooth,thresh)
    
    if strcmp(thresh.units,'absolute')
        threshInt = thresh.level;
    elseif strcmp(thresh.units,'SD')
        threshInt = thresh.level*thresh.sd;
    end
    
    stack2 = smooth - threshInt > 0;
    npix_above_thresh = sum(sum(sum(stack2)));
    
    clear('threshInt');
end