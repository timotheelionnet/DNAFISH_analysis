function [new_ctr,ROIlimits] = compute_ROI_boundariesYi(stack,spot_ctr,cutwidth,thickness,ROIsize)

%check input consistency
if ndims(stack) == 3 && numel(spot_ctr) < 3
    disp('cannot compute ROI boundaries, incorrect ROI center coordinates');
end

%compute half length of ROI square/cube
if strcmp(ROIsize,'small')
    halflength = cutwidth(1);
    if ndims(stack) == 3
        halflength(2) = cutwidth(2);
    end
elseif strcmp(ROIsize,'large')
    halflength = cutwidth(1)+thickness;
    if ndims(stack) == 3
        halflength(2) = cutwidth(2)+thickness;
    end
end

%input size
if ndims(stack) == 3
    [nx,ny,nz]= size(stack);
    xc = spot_ctr(1); yc = spot_ctr(2); zc = spot_ctr(3);
else
    [nx,ny]= size(stack);
    xc = spot_ctr(1); yc = spot_ctr(2); 
end

%compute ROI limits
xmin = ceil(xc) - halflength(1);
xmin = max(1,xmin);
ymin = ceil(yc) - halflength(1);
ymin = max(1,ymin);

xmax = ceil(xc) + halflength(1);
xmax = min(nx,xmax);
ymax = ceil(yc) + halflength(1);
ymax = min(ny,ymax);

ROIlimits = [xmin,ymin;xmax,ymax];

if ndims(stack) == 3
    zmin = ceil(zc) - halflength(2);
    zmin = max(1,zmin);
    
    zmax = ceil(zc) + halflength(2);
    zmax = min(nz,zmax);
    
    ROIlimits = [ROIlimits,[zmin;zmax]];
end

%compute coordinates of spot center in new region
new_ctr(1) = spot_ctr(1) - xmin + 1;
new_ctr(2) = spot_ctr(2) - ymin + 1;
if ndims(stack) == 3
    new_ctr(3) = spot_ctr(3) - zmin + 1;
end

end