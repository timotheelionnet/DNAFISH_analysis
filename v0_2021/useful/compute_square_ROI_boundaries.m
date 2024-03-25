function [new_ctr,ROIlimits] = compute_square_ROI_boundaries(imgsize,spot_ctr,square_size_in_pix)
%INPUT
%imgsize: size of image [nx,ny] or [ny,ny,nz]
%spot_ctr: position of spot center in pixels: [xc,yc] or [yc,yc,zc]
%square_size_in_pix: size of the square or cube around each spot
%[2*dx+1,2*dy+1] or [2*dx+1,2*dy+1, 2*dz+1] 

%OUTPUT
%new_ctr are the coordinates of spot_ctr relative to the new ROI
%ROIlimits [xmin,ymin;xmax,ymax;] or [xmin,ymin,zmin;xmax,ymax,zmax;] the
%indices of the corner voxels of the ROI

%check input consistency
if numel(imgsize) ~= numel(spot_ctr) || numel(imgsize) ~= numel(square_size_in_pix)
    dispwin('error','cannot compute ROI boundaries, mismatch between ROI/center/dimension coordinates');
end

%compute half length of ROI square/cube
halflength = (square_size_in_pix-1)/2;


%input size
if  numel(imgsize) ~= 2 && numel(imgsize) ~= 3
    disp('image size not recognized');
    new_ctr = [];
    ROIlimits = [];
    return
    
elseif numel(imgsize) == 3
    [nx,ny,nz] = deal(imgsize(1),imgsize(2),imgsize(3));
    xc = spot_ctr(1); 
    yc = spot_ctr(2); 
    zc = spot_ctr(3);
    
elseif numel(imgsize) == 2
    [nx,ny] = deal(imgsize(1),imgsize(2));
    xc = spot_ctr(1); 
    yc = spot_ctr(2); 
end

%compute ROI limits
xmin = ceil(xc) - halflength(1);
xmin = max(1,xmin);
ymin = ceil(yc) - halflength(2);
ymin = max(1,ymin);

xmax = ceil(xc) + halflength(1);
xmax = min(nx,xmax);
ymax = ceil(yc) + halflength(2);
ymax = min(ny,ymax);

ROIlimits = [xmin,ymin;xmax,ymax];

if numel(imgsize) == 3
    zmin = round(zc) - halflength(3);
    zmin = max(1,zmin);
    
    zmax = round(zc) + halflength(3);
    zmax = min(nz,zmax);
    
    ROIlimits = [ROIlimits,[zmin;zmax]];
end

%compute coordinates of spot center in new region
new_ctr(1) = spot_ctr(1) - xmin + 1;
new_ctr(2) = spot_ctr(2) - ymin + 1;
if numel(imgsize) == 3
    new_ctr(3) = spot_ctr(3) - zmin + 1;
end

end

