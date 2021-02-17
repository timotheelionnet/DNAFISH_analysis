function plot_current_z_stack_3colors(fh,ha,hImin,hImax,stack1,img2,img3,z,hoverlay)
%stack1: green
%img2: red
%img3: blue
%% parsing dimensions of input stacks    
    if ndims(stack1)==3
        [nx ny nz] = size(stack1);
        if z<=0 || z >nz
            disp('z out of bounds in plot_current_z_stack!');
        end
        stack1 = stack1(:,:,z); 
    
    elseif ndims(stack1) == 2
        [nx ny] = size(stack1);
    else
        disp('wrong stack size');
    end
    
    if ndims(img2) == 3
        img2 = img2(:,:,z);
    elseif img2 == 0
        img2 = zeros(nx,ny);
    end
    if ndims(img3) == 3
        img3 = img3(:,:,z);
    elseif img3 == 0
        img3 = zeros(nx,ny);    
    end
    
    Imin = str2num(get(hImin,'String'));
    Imax = str2num(get(hImax,'String'));
    
    cur_slice = zeros(nx,ny,3);

%% filling channels  

    %red channel: overlay or nothing
    is_overlay_on = 1 - get(hoverlay,'Value');
    if is_overlay_on
        cur_slice(1:nx,1:ny,1) = img2(1:nx,1:ny);
    else
        cur_slice(1:nx,1:ny,1) = 0;
    end
    
    %green channel: data
    if(Imin >= Imax)
        cur_slice(1:nx,1:ny,2) = 0;
    else        
        cur_slice(1:nx,1:ny,2) = min(max( (double(stack1(1:nx,1:ny))-Imin) / (Imax-Imin),0),1); 
    end

    %blue channel: overlay or nothing
    is_overlay_on = 1 - get(hoverlay,'Value');
    if is_overlay_on
        cur_slice(1:nx,1:ny,3) = img3(1:nx,1:ny);
    else
        cur_slice(1:nx,1:ny,3) = 0;
    end

%%    
    set(fh,'CurrentAxes',ha);
    imagesc(cur_slice); 
    xlim(ha,[1,max(nx,ny)]);
    ylim(ha,[1,max(nx,ny)]);

clear('cur_slice','nx','ny','nz','Imin','Imax','is_overlay_on','newXLim','newYLim','stack1','stack2');
clear('fh','ha','hImin','hImax','stack1','stack2','z','hoverlay');
end