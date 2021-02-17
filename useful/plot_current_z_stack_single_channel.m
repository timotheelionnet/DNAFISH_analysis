function plot_current_z_stack_single_channel(fh,ha,hImin,hImax,stack1,z)
    
%% parsing dimensions of input stack    
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
    
    Imin = str2num(get(hImin,'String'));
    Imax = str2num(get(hImax,'String'));

    cur_slice = zeros(nx,ny);

%% filling data 

    %green channel: data
    if(Imin >= Imax)
        cur_slice(1:nx,1:ny) = 0;
    else
        cur_slice(1:nx,1:ny) = min( max((stack1(1:nx,1:ny)-Imin)/(Imax-Imin),0) , 1);
    end

%%    
    set(fh,'CurrentAxes',ha);
    
    imagesc(cur_slice); colormap(gray);
    xlim(ha,[1,max(nx,ny)]);
    ylim(ha,[1,max(nx,ny)]);
    
    %set(ha,'XLim',newXLim);
    %set(ha,'YLim',newYLim);

clear('cur_slice','nx','ny','nz','Imin','Imax','newXLim','newYLim','stack1');
clear('fh','ha','hImin','hImax','stack1','z');