function plot_zoom_single_channel(x,y,z,width,zoomh,stack1,fh,hImin,hImax)
    
    set(fh,'CurrentAxes',zoomh);
    if ndims(stack1) == 3
        [nx ny nz] = size(stack1);
        if z<=0 || z >nz
            disp('z out of bounds in plot_current_z_stack!');
        end
        stack1 = stack1(:,:,z); 
        
    elseif ndims(stack1) ==2
        [nx ny] = size(stack1);
    end
   
    Imin = str2num(get(hImin,'String'));
    Imax = str2num(get(hImax,'String'));

    %computing the value of the local slice 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cur_slice = zeros(nx,ny);

    if(Imin >= Imax)
        cur_slice(1:nx,1:ny) = 0;
    else
        cur_slice(1:nx,1:ny) = min(1, max( (stack1(1:nx,1:ny)-Imin) / (Imax-Imin),0));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %chosing the window corresponding to the cursor position
    if 2*width+1>nx
        if 2*width+1>ny
            imagesc(cur_slice);
        elseif (y - width >= 1) && (y + width <= ny)
            imagesc([y-width y+width],[1 nx],cur_slice(1:nx,y-width:y+width));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            imagesc([yy-width,yy+width],[1,nx],cur_slice(1:nx,yy-width:yy+width));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            imagesc([yy-width,yy+width],[1,nx],cur_slice(1:nx,yy-width:yy+width));
        end
    elseif (x - width >= 1) && (x + width <= nx)
        if 2*width+1>ny
            imagesc([1,ny],[x-width,x+width],cur_slice(x-width:x+width,1:ny));
        elseif (y - width >= 1) && (y + width <= ny)
            imagesc([y-width,y+width],[x-width,x+width],cur_slice(x-width:x+width,y-width:y+width));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            imagesc([yy-width,yy+width],[x-width,x+width],cur_slice(x-width:x+width,yy-width:yy+width));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            imagesc([yy-width,yy+width],[x-width,x+width],cur_slice(x-width:x+width,yy-width:yy+width));
        end
    elseif (x - width < 1) && (x + width <= nx) 
        xx = width+1;
        if 2*width+1>ny
            imagesc([1,ny],[xx-width,xx+width],cur_slice(xx-width:xx+width,1:ny));
        elseif (y - width >= 1) && (y + width <= ny)
            imagesc([y-width,y+width],[xx-width,xx+width],cur_slice(xx-width:xx+width,y-width:y+width));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            imagesc([yy-width,yy+width],[xx-width,xx+width],cur_slice(xx-width:xx+width,yy-width:yy+width));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            imagesc([yy-width,yy+width],[xx-width,xx+width],cur_slice(xx-width:xx+width,yy-width:yy+width));
        end
    elseif (x + width > nx) && (x - width >= 1)
        xx = nx - width;
        if 2*width+1>ny
            imagesc([1,ny],[xx-width,xx+width],cur_slice(xx-width:xx+width,1:ny));
        elseif (y - width >= 1) && (y + width <= ny)
            imagesc([y-width,y+width],[xx-width,xx+width],cur_slice(xx-width:xx+width,y-width:y+width));
        elseif (y - width < 1) && (y + width <= ny)
            yy = width+1;
            imagesc([yy-width,yy+width],[xx-width,xx+width],cur_slice(xx-width:xx+width,yy-width:yy+width));
        elseif (y + width > ny) && (y - width >= 1)
            yy = ny - width;
            imagesc([yy-width,yy+width],[xx-width,xx+width],cur_slice(xx-width:xx+width,yy-width:yy+width));
        end
    end 
    colormap(gray);
    clear('cur_slice','nx','ny','nz','Imin','Imax','stack1');
    clear('xx','yy');
    clear('x','y','z','width','zoomh','stack1','fh','hImin','hImax');
end