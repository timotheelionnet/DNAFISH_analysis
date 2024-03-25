function p = smooth_image_and_subtract_background5(p)
    %p: parameters structure
    
    
    if p.numdim == 3
        disp('3D smoothing of the image ...');
        p.smooth = bpass_filter_3D_fast_memory_safe3(p.data,p.filter.nlo,p.sigma_xy,p.sigma_z,p.filter.nhi,p.filter.width,p.filter.numdim);
        
    elseif p.numdim == 2
        disp('2D smoothing of the image ...');
        p.smooth = bpass_filter_2D_fourier(p.data,p.filter.nlo,p.sigma_xy,p.filter.nhi,p.filter.width); 
    else
        disp('smoothing only possible on 2d or 3d stacks');
        p.smooth = 0;
        return;
    end
    
    %I subtract the background
    x= mean(p.smooth(:));
    p.smooth = p.smooth - x;
    if p.numdim == 3
        p.thresh.sd = std_stack(p.smooth);
    else
        p.thresh.sd = std_im(p.smooth);
    end
    clear('x');
end

