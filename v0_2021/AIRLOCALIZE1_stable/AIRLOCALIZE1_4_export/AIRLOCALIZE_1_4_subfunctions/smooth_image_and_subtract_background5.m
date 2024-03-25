function p = smooth_image_and_subtract_background5(p)
    %p: parameters structure
    %version 5 of the function uses only the parameters
    %p.data: array to filter
    %p.filter.nlo:
    %p.filter.nhi
    %p.sigma_xy
    
    %the follwoing parameters are not in use anymore:
    %p.filter
    %p.sigma_z
    %p.filter_width
    %p.filter.numdim
    
    %note that the default values for the cutoffs are different for the DOG
    %filter compared to the previous Fourier based filter
    
    if p.numdim == 3 || p.numdim == 2
        fh = dispwin('progress','smoothing image ...');
        %factors 1.5 is a heuristic attempt to reproduce results
        %from the Fourier filter given the same parameters
        p.smooth = 1.5*DOGfilter(p.data,p.filter.nhi,p.filter.nlo*p.sigma_xy,[],[]);
    else
        fh = dispwin('data type error','smoothing only possible on 2d images or 3d stacks');
        p.smooth = 0;
        close(fh);
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
    close(fh);
end

