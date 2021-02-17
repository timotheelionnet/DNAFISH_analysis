function im2 = adjust_size_stack(im,nx,ny)

if ndims(im) == 2
    [mx,my] = size(im);
    mz = 1;
elseif ndims(im) == 3
    [mx,my,mz] = size(im);
else
    disp('images should be 2d or 3d');
end

if mx>= nx && my >= ny
    im2 = im;
    return
end

im2 = zeros(max(mx,nx),max(my,ny),mz);

if mx>= nx && my < ny
    dy = ny - my;
    dy1 = round(dy/2);
    
    im2(1:mx,dy1:dy1+my,1:mz) = im;   
end

if my>= ny && mx < nx
    dx = nx - mx;
    dx1 = round(dx/2);
    
    im2(dx1:dx1+mx,1:my,1:mz) = im;   
end

if my < ny && mx < nx
    dx = nx - mx;
    dx1 = round(dx/2);
    dy = ny - my;
    dy1 = round(dy/2);
    
    im2(dx1:dx1+mx,dy1:dy1+my,1:mz) = im;   
end



end