function [pl,stack_bg_corr,new_ctr,ROIlimits] = subtract_median(stack,spot_ctr,cutwidth,thickness,ROIsize,median_type)

[new_ctr,ROIlimits] = compute_ROI_boundaries(stack,spot_ctr,cutwidth,thickness,ROIsize);

if ndims(stack) == 3
    stack_bg_corr = stack(...
        ROIlimits(1,1):ROIlimits(2,1),...
        ROIlimits(1,2):ROIlimits(2,2),...
        ROIlimits(1,3):ROIlimits(2,3));
else
    stack_bg_corr = stack(...
        ROIlimits(1,1):ROIlimits(2,1),...
        ROIlimits(1,2):ROIlimits(2,2));
end

if strcmp(median_type,'local')
    stack_bg_corr = stack_bg_corr - median(double(stack_bg_corr(:)));
elseif strcmp(median_type,'global')
    stack_bg_corr = stack_bg_corr - median(double(stack(:)));
end

pl = [];
end