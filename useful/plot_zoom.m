function plot_zoom(x,y,z,width,zoomh,stack1,stack2,fh,hImin,hImax,hoverlay)
    
    set(fh,'CurrentAxes',zoomh);
    if ndims(stack1) == 3
        [nx ny nz] = size(stack1);
        if z<=0 || z >nz
            disp('z out of bounds in plot_current_z_stack!');
        end
        stack1 = stack1(:,:,z); stack2 = stack2(:,:,z);
        
    elseif ndims(stack1) ==2
        [nx ny] = size(stack1);
    end
   
    Imin = str2num(get(hImin,'String'));
    Imax = str2num(get(hImax,'String'));

    %computing the value of the local slice 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cur_slice = zeros(nx,ny,3);

    %red channel: overlay or nothing
    is_overlay_on = 1 - get(hoverlay,'Value');
    if is_overlay_on
        cur_slice(1:nx,1:ny,1) = stack2(1:nx,1:ny);
    else
        cur_slice(1:nx,1:ny,1) = 0;
    end

    %green channel: data
    if(Imin >= Imax)
        cur_slice(1:nx,1:ny,2) = 0;
    else
        cur_slice(1:nx,1:ny,2) = min(max( (double(stack1(1:nx,1:ny))-Imin) / (Imax-Imin),0),1);
    end

    %nothing in the blue channel
    cur_slice(1:nx,1:ny,3) = 0; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %chosing the window corresponding to the cursor position
    if 2*width+1>nx
        if 2*width+1>ny
            image(cur_slice);
        elseif (y - width >= 1) && (y + width <= ny)
            image([y-width y+width],[1 nx],cur_slice(1:nx,y-width:y+width,1:3));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            image([yy-width,yy+width],[1,nx],cur_slice(1:nx,yy-width:yy+width,1:3));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            image([yy-width,yy+width],[1,nx],cur_slice(1:nx,yy-width:yy+width,1:3));
        end
    elseif (x - width >= 1) && (x + width <= nx)
        if 2*width+1>ny
            image([1,ny],[x-width,x+width],cur_slice(x-width:x+width,1:ny,1:3));
        elseif (y - width >= 1) && (y + width <= ny)
            image([y-width,y+width],[x-width,x+width],cur_slice(x-width:x+width,y-width:y+width,1:3));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            image([yy-width,yy+width],[x-width,x+width],cur_slice(x-width:x+width,yy-width:yy+width,1:3));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            image([yy-width,yy+width],[x-width,x+width],cur_slice(x-width:x+width,yy-width:yy+width,1:3));
        end
    elseif (x - width < 1) && (x + width <= nx) 
        xx = width+1;
        if 2*width+1>ny
            image([1,ny],[xx-width,xx+width],cur_slice(xx-width:xx+width,1:ny,1:3));
        elseif (y - width >= 1) && (y + width <= ny)
            image([y-width,y+width],[xx-width,xx+width],cur_slice(xx-width:xx+width,y-width:y+width,1:3));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            image([yy-width,yy+width],[xx-width,xx+width],cur_slice(xx-width:xx+width,yy-width:yy+width,1:3));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            image([yy-width,yy+width],[xx-width,xx+width],cur_slice(xx-width:xx+width,yy-width:yy+width,1:3));
        end
    elseif (x + width > nx) && (x - width >= 1)
        xx = nx - width;
        if 2*width+1>ny
            image([1,ny],[xx-width,xx+width],cur_slice(xx-width:xx+width,1:ny,1:3));
        elseif (y - width >= 1) && (y + width <= ny)
            imagesc([y-width,y+width],[xx-width,xx+width],cur_slice(xx-width:xx+width,y-width:y+width,1:3));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            imagesc([yy-width,yy+width],[xx-width,xx+width],cur_slice(xx-width:xx+width,yy-width:yy+width,1:3));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            imagesc([yy-width,yy+width],[xx-width,xx+width],cur_slice(xx-width:xx+width,yy-width:yy+width,1:3));
        end
    end 
    
    clear('cur_slice','nx','ny','nz','Imin','Imax','is_overlay_on','stack1','stack2');
    clear('xx','yy');
    clear('x','y','z','width','zoomh','stack1','stack2','fh','hImin','hImax','hoverlay');
end