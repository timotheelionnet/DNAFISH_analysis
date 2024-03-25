function [x0,y0,z0,N0,err0,dist,it] = apply_gaussian_mask_fixed_position(Im,spot_ctr,p)
% runs a 2D/3D gaussian mask algorithm to localize and quantify the intensity
% of a fluorescent spot
%Im is the image/stack
%all units in pix with origin at the edge of pixel i.e. leftmost corner
%center corresponds to y = 0.5
%spot_ctr = [x y z] : the guess center coordinates (or [x,y] in 2D)
%p is the structure holding the fit parameters
%p.type : 'gaussian' or 'integrated gaussian' or 'integrated gaussian std
%z' (last option only in 3D)
%p.sigma_xy : PSF width (in pixels)
%p.sigma_z : PSF height (in pixels - ignored in 2D)
%p.cutsize : range (in PSF width units) of the region around the spot used for fitting. 
%p.maxcount is the maximum number of iterations of the equations allowed
%before convergence
%p.tol is the tolerance on the convergence (in lateral pixel dimension units)

%% parse input
Im = double(Im);
xs = double(spot_ctr(1)); 
ys = double(spot_ctr(2)); 
sxy = double(p.sigma_xy); 
tol = double(p.tol);
maxcount = p.maxcount;
cutsize = double(p.cutsize);
cutwidth = cutsize*sxy;
dx = 1;
dy = 1;
if ndims(Im) == 3
    [nx,ny,nz] = size(Im);
    zs = double(spot_ctr(3));
    sz = double(p.sigma_z);
    cutwidth(2) = cutsize*sz;
    dz = 1;
else
    [nx,ny] = size(Im);
    z0 = 0;
end


%% compute boundaries of the ROI over which the mask is applied
if ndims(Im) == 3
    [~,ROIlimits] = compute_ROI_boundaries(Im,[xs,ys,zs],cutwidth,0,'small');
    %voxel indices (integer)
    xp_min = ceil(ROIlimits(1,1)); 
    xp_max = ceil(ROIlimits(2,1));
    yp_min = ceil(ROIlimits(1,2)); 
    yp_max = ceil(ROIlimits(2,2));
    zp_min = round(ROIlimits(1,3)); 
    zp_max = round(ROIlimits(2,3));
    [xp,yp,zp] = ndgrid(xp_min:xp_max,yp_min:yp_max,zp_min:zp_max);
    sIm = double(Im(xp_min:xp_max,yp_min:yp_max,zp_min:zp_max));
else
    [~,ROIlimits] = compute_ROI_boundaries(Im,[xs,ys],cutwidth,0,'small');
    %pixel indices (integer)
    xp_min = ceil(ROIlimits(1,1)); 
    xp_max = ceil(ROIlimits(2,1));
    yp_min = ceil(ROIlimits(1,2)); 
    yp_max = ceil(ROIlimits(2,2));
    [xp,yp] = ndgrid(xp_min:xp_max,yp_min:yp_max);
    sIm = double(Im(xp_min:xp_max,yp_min:yp_max));
end

%% apply mask

    if ndims(Im) == 3
        switch p.type
            case 'gaussian'
                intensity = intensity_gaussian3D(xp,yp,zp,xs,ys,zs,sxy,sz,dx,dy,dz);               
            case 'integrated gaussian' 
                intensity = intensity_integrated_gaussian3D(xp,yp,zp,xs,ys,zs,sxy,sz,dx,dy,dz);
            case 'integrated gaussian std z' 
                intensity = intensity_integrated_gaussian3D_stdZ(xp,yp,zp,xs,ys,zs,sxy,sz,dx,dy,dz);                
        end
        intsum = intensity .* sIm;
        intsum = sum(intsum(:));
        
        sumsum = intensity.*intensity;
        sumsum = sum(sumsum(:));
        
        sumx = (double(xp-0.5)).*intensity.*sIm;
        sumx = sum(sumx(:));
        
        sumy = (double(yp-0.5)).*intensity.*sIm;
        sumy = sum(sumy(:));
        
        sumz = double(zp).*intensity.*sIm;
        sumz = sum(sumz(:));
        
        if intsum <= 0 || sumsum == 0
            [x0,y0,z0,N0] = deal(0,0,0,0);
        else
            
            x0 = double(sumx) / double(intsum);
            y0 = double(sumy) / double(intsum);
            z0 = double(sumz) / double(intsum);
            N0 = double(intsum) / double(sumsum);  
%             if   location_out_of_ROI([x0,y0,z0],[xp_min,xp_max,yp_min,yp_max,zp_min,zp_max],[dx,dy,dz])
%                    [x0,y0,z0,N0] = deal(0,0,0,0);
%             end 
        end
  
    elseif ndims(Im) == 2
        switch p.type
            case 'gaussian'
                intensity = intensity_gaussian2D(xp,yp,xs,ys,sxy,dx,dy);
            case 'integrated gaussian' 
                intensity = intensity_integrated_gaussian2D(xp,yp,xs,ys,sxy,dx,dy);
            case 'integrated gaussian std z' 
                intensity = intensity_integrated_gaussian2D(xp,yp,xs,ys,sxy,dx,dy);
        end
        intsum = intensity .* sIm;
        intsum = sum(intsum(:));
        
        sumsum = intensity.*intensity;
        sumsum = sum(sumsum(:)) ;
        
        sumx = (double(xp-0.5)).*intensity.*sIm;
        sumx = sum(sumx(:));
        
        sumy = (double(yp-0.5)).*intensity.*sIm;
        sumy = sum(sumy(:));
        
        if intsum <= 0 || sumsum == 0
            [x0,y0,N0] = deal(0,0,0);
        else
            x0 = double(sumx) / double(intsum);
            y0 = double(sumy) / double(intsum);
            N0 = double(intsum) / double(sumsum); 
%             if location_out_of_ROI([x0,y0],[xp_min,xp_max,yp_min,yp_max],[dx,dy])
%                        [x0,y0,N0] = deal(0,0,0);
%             end      
        end
    end      


%% error computation
if ndims(Im) == 3
    err0 = N0*double(intensity) - double(Im(xp_min:xp_max,yp_min:yp_max,zp_min:zp_max));
    err0 = err0.^2;
    err0 = sqrt(sum(err0(:)));    
elseif ndims(Im) == 2
    err0 = N0*double(intensity) - double(Im(xp_min:xp_max,yp_min:yp_max));
    err0 = err0.^2;
    err0 = sqrt(sum(err0(:)));
end

%% convert PSF prefactor into integrated intensity
if ndims(Im) == 3
    switch p.type
        case 'integrated gaussian'
           N0 = 8 * N0; 
        case 'integrated gaussian std z'
           x = (floor(x0 - 3*sxy): ceil(x0 + 3*sxy )  );
           y = (floor(y0 - 3*sxy): ceil(y0 + 3*sxy ) );
           z = (floor(z0 - 3*sz): ceil(z0 + 3*sz )  );
           [yy,xx,zz] = meshgrid(y,x,z);
           xx = ceil(xx)-0.5; 
           yy = ceil(yy)-0.5;
           zz = round(zz);
           Itot = intensity_integrated_gaussian3D_stdZ(xx,yy,zz,xs,ys,zs,sxy,sz,1,1,1);
           N0 = sum(Itot(:))*N0; 
        case 'gaussian'
           x = (floor(x0 - 3*sxy): ceil(x0 + 3*sxy )  );
           y = (floor(y0 - 3*sxy): ceil(y0 + 3*sxy )  );
           z = (floor(z0 - 3*sz): ceil(z0 + 3*sz )  );
           [yy,xx,zz] = meshgrid(y,x,z);
           xx = ceil(xx)-0.5; 
           yy = ceil(yy)-0.5;
           zz = round(zz);
           Itot = intensity_gaussian3D(xx,yy,zz,xs,ys,zs,sxy,sz,1,1,1);
           N0 = sum(Itot(:))*N0;       
    end
else
    switch p.type
        case 'integrated gaussian'
           N0 = 4 * N0; 
        case 'integrated gaussian std z'
           N0 = 4 * N0;  
        case 'gaussian'
           x = (floor(x0 - 3*sxy): ceil(x0 + 3*sxy )  );
           y = ( floor(y0 - 3*sxy): ceil(y0 + 3*sxy )  );
           [yy,xx] = meshgrid(y,x);
           xx = ceil(xx)-0.5; 
           yy = ceil(yy)-0.5;
           Itot = intensity_gaussian2D(xx,yy,xs,ys,sxy,1,1);
           N0 = sum(Itot(:))*N0;    
    end
end

end

function res = location_out_of_ROI(location,boundaries,vox_size)
    %test if the location found by the mask is within the ROI
    %boundaries in pizel units
    %x0 y0 z0 refer to positions with origin on the edge of the pixel.
 
    if numel(location) == 3
        x0 = location(1);
        y0 = location(2);
        z0 = location(3);
        xp_min = boundaries(1);
        xp_max = boundaries(2);
        yp_min = boundaries(3);
        yp_max = boundaries(4);
        zp_min = boundaries(5);
        zp_max = boundaries(6);
        dx = vox_size(1);
        dy = vox_size(2);
        dz = vox_size(3);
    else
        x0 = location(1);
        y0 = location(2);
        xp_min = boundaries(1);
        xp_max = boundaries(2);
        yp_min = boundaries(3);
        yp_max = boundaries(4);
        dx = vox_size(1);
        dy = vox_size(2);
    end
    
    res = 0;    
    if x0/dx <xp_min-1
        res = 1;
        return
    end
    if x0/dx >xp_max
        res = 1;
        return
    end
    
    if y0/dy <yp_min -1
        res = 1;
        return
    end
    if y0/dy >yp_max
        res = 1;
        return
    end
    
    if numel(location) == 3
        if z0/dz <zp_min -1
            res = 1;
            return
        end
        if z0/dz >zp_max
            res = 1;
            return
        end
    end
    
end