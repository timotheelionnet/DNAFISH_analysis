function [xp_min yp_min zp_min xp_max yp_max zp_max] = get_ROI_boundaries_in_pix(xs,ys,zs,cutsize,bgcutsize,sxy,sz,dx,dy,dz,nx,ny,nz)
%get the boundaries of an ROI centered on a point  of coordinates [xs ys zs] in nm
%cutsize is in units of the PSF width (resp sxy and sz)
%includes a bgcutsize for background correction in units of pixels
%dx dy dz are voxel dimensions.
display_message = 0;

xmin = xs - cutsize*sxy - bgcutsize*dx;
xmax = xs + cutsize*sxy + bgcutsize*dx;
ymin = ys - cutsize*sxy - bgcutsize*dy;
ymax = ys + cutsize*sxy + bgcutsize*dy;
zmin = zs - cutsize*sz  - bgcutsize*dz;
zmax = zs + cutsize*sz  + bgcutsize*dz;
            
%same locations in pixel units
xp_min = max( ceil(xmin/dx), 1);
if xp_min == 1 && display_message==1
    disp('xmin out of bounds, careful');
end
xp_max = min( ceil(xmax/dx), nx );
if xp_max == nx  && display_message==1
    disp('xmax out of bounds, careful');
end
yp_min = max( ceil(ymin/dy), 1);
if yp_min == 1 && display_message==1
    disp('ymin out of bounds, careful');
end
yp_max = min( ceil(ymax/dy), ny );
if yp_max == ny && display_message==1
    disp('ymax out of bounds, careful');
end
zp_min = max( ceil(zmin/dz), 1);
if zp_min == 1 && display_message==1
    disp('zmin out of bounds, careful');
end
zp_max = min( ceil(zmax/dz), nz );
if zp_max == nz && display_message==1
    disp('zmax out of bounds, careful');
end
