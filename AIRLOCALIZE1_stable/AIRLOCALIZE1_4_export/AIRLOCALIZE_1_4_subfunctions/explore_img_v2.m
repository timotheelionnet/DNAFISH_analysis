function varargout = explore_img_v2(varargin)
% VIEW_STACK 

%one argument: one color channel (data)
%two arguments: the first one is the data (green channel)
%the second one is the overlay (red channel) - typically the position of
%the spots located by some detection algorithm
%three arguments: the first one is the data (green channel)
%the second one is the overlay (red channel) - typically the position of
%the spots located by some detection algorithm
%the third one is a string specifying the shape of the object locating the spots in the overlay channel


%% parsing arguments 
%checking that we have at least one 3D stack or load one. this is the data (green channel)
if nargin == 0
    [fname,sourcedir,fidx] = uigetfile('*.tif;*.stk;*.lsm','Select Image File to Estimate PSF');
    if fidx == 0
        varargout(1) = 'cancel'; 
        clear('fname','sourcedir','fidx');
        return; 
    end
    stack1 = timtiffread(fullfile(sourcedir,fname));
    stack1name = fname;
    
elseif nargin >= 1
    stack1 = varargin{1};
    stack1 = double(stack1); 
    stack1name = inputname(1);
    if nargin >= 2
        if isstruct(varargin{2})
            p = varargin{2};
        else
            p = 0;
        end
    else
        p = 0;    
    end
end

if ndims(stack1) ~= 2 || size(stack1,2)==1
    fprintf('data should be 2 dimensional.\n');
    return;
end

%% main figure
fh = figure(...
              'Units','characters',...
              'MenuBar','none',...
              'Toolbar','none',...
              'NumberTitle','off',...
              'Position',[20 3 200 55],...
              'Visible','off');  


[nx ny] = size(stack1);
setappdata(fh,'nx',nx);
setappdata(fh,'ny',ny);
xfreeze = round(nx/2); 
yfreeze = round(ny/2); 
setappdata(fh,'xfreeze',xfreeze);
setappdata(fh,'yfreeze',yfreeze);
setappdata(fh,'fitres',[0,0,0,0,0]);
fh = set_figure_and_panel_for_explore_img(stack1,stack1name,fh);
handles = guihandles(fh);

%%%%%%%%%%%%%%%%%%%%%%%%%%% set all the callback functions %%%%%%%%%%%%%%%%%            
set(fh,'WindowButtonMotionFcn', {@viewlocaldata,stack1,fh});
set(fh,'WindowButtonUpFcn', {@change_state,fh});   
set(handles.hclose,'Callback', {@clean_close,fh});   
set(handles.hImin,'Callback',{@change_Iminmax,stack1,fh});
set(handles.hImax,'Callback',{@change_Iminmax,stack1,fh});
set(handles.hglob_auto,'Callback',{@global_auto_adjust,stack1,fh});
set(handles.hgauss_fit,'Callback',{@local_2D_gaussian_fit_from_viewer,stack1,p,fh});
set(handles.hrecord,'Callback',{@record_fit_results,stack1,fh});
set(handles.hxplot_range,'Callback',{@change_local_range,stack1,fh});
set(handles.hyplot_range,'Callback',{@change_local_range,stack1,fh});

set(fh,'SelectionType','alt');
set(fh,'CurrentAxes',handles.ha);
    
plot_current_z_stack_single_channel(fh,stack1,1);
set(fh,'Visible','on');

uiwait;

varargout{1} = getappdata(fh,'fitres');
close(fh);
drawnow;
%% housekeeping
clear('nx','ny','nz','stack1','stack1name','xfreeze','yfreeze');
clear('fh','ha','handles');

end

%%%%%%%%%%%%%%%%%%%%%%%%%% Callback functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_current_z_stack_single_channel(fh,stack1,z)

handles = guihandles(fh);    
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
    
    Imin = str2num(get(handles.hImin,'String'));
    Imax = str2num(get(handles.hImax,'String'));

    cur_slice = zeros(nx,ny);

%% filling data 

    %green channel: data
    if(Imin >= Imax)
        cur_slice(1:nx,1:ny) = 0;
    else
        cur_slice(1:nx,1:ny) = min( max((double(stack1(1:nx,1:ny))-Imin)/(Imax-Imin),0) , 1);
    end

%%    
    ha = handles.ha;
    set(fh,'CurrentAxes',handles.ha);
    
    %ha = imagesc(cur_slice,'XLim',[1,max(nx,ny)],'YLim',[1,max(nx,ny)]); colormap(gray);
    imagesc(cur_slice); colormap(gray);
    xlim([1,max(nx,ny)]);
    ylim([1,max(nx,ny)]);
    set(ha,'Tag','ha');
    
    %set(ha,'XLim',[1,max(nx,ny)] );
    %set(ha,'YLim',[1,max(nx,ny)] );

clear('cur_slice','nx','ny','nz','Imin','Imax','newXLim','newYLim','stack1');
clear('stack1','handles','ha');
end

function plot_zoom_single_channel(x,y,z,width,stack1,fh)
    
    handles = guihandles(fh); 
    zoomh = handles.zoomh;
    set(fh,'CurrentAxes',handles.zoomh);
    if ndims(stack1) == 3
        [nx ny nz] = size(stack1);
        if z<=0 || z >nz
            disp('z out of bounds in plot_current_z_stack!');
        end
        stack1 = stack1(:,:,z); 
        
    elseif ndims(stack1) ==2
        [nx ny] = size(stack1);
    end
   
    Imin = str2num(get(handles.hImin,'String'));
    Imax = str2num(get(handles.hImax,'String'));

    %computing the value of the local slice 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cur_slice = zeros(nx,ny);

    if(Imin >= Imax)
        cur_slice(1:nx,1:ny) = 0;
    else
        cur_slice(1:nx,1:ny) = min(1, max( (double(stack1(1:nx,1:ny))-Imin) / (Imax-Imin),0));
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
    set(zoomh,'Tag','zoomh');
    colormap(gray);
    
    clear('cur_slice','nx','ny','nz','Imin','Imax','stack1','handles');
    clear('xx','yy');
    clear('x','y','z','width','zoomh','stack1','fh','hImin','hImax');
end

function change_local_range(src,eventdata,stack1,fh)
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    plot_local_data_single_channel(x,y,1,fh,stack1);
    clear('x','y','z','handles');
end

function record_fit_results(src,eventdata,stack1,fh)
    handles = guihandles(fh);
    fit_res = getappdata(fh,'fitres');
    cur_res(1,1) = str2double(get(handles.hx0,'String'));
    cur_res(1,2) = str2double(get(handles.hy0,'String'));
    cur_res(1,3) = str2double(get(handles.hI0,'String'));
    cur_res(1,4) = str2double(get(handles.hbg0,'String'));
    cur_res(1,5) = str2double(get(handles.hsxy0,'String'));
    
    n = str2double(get(handles.hFitsRecorded,'String'));
    if n == 0
        setappdata(fh,'fitres',cur_res);
    else
        setappdata(fh,'fitres',[fit_res;cur_res]);
    end
    set(handles.hFitsRecorded,'String',num2str(n+1));
    
    clear('handles','fit_res','cur_res','n');
end

function clean_close(src,eventdata,fh)
    uiresume;
end

function local_2D_gaussian_fit_from_viewer(src,eventdata,stack1,p,fh)

    xfreeze = getappdata(fh,'xfreeze');
    yfreeze = getappdata(fh,'yfreeze');
    handles = guihandles(fh);
     
    Gaussout = gaussian_fit_local(stack1,[xfreeze,yfreeze],p,1);
    
    set(handles.hI0,'String',num2str(Gaussout(1)));
    set(handles.hbg0,'String',num2str(Gaussout(2)));
    set(handles.hsxy0,'String',num2str(Gaussout(3)));
    set(handles.hx0,'String',num2str(Gaussout(4)));
    set(handles.hy0,'String',num2str(Gaussout(5)));
    
    plot_local_data_and_fit_single_channel_bgcorr(xfreeze,yfreeze,1,fh,stack1,p,Gaussout);
    clear('xfreeze','yfreeze','handles','z','p','Gaussout','resnorm');
end

function change_state(src,eventdata,fh)
    %switches between the snap and grab mode
    %hstate value = 1: grab
    %hstate value = 2: snap
    
    xfreeze = getappdata(fh,'xfreeze');
    yfreeze = getappdata(fh,'yfreeze');
    nx = getappdata(fh,'nx');
    ny = getappdata(fh,'ny');
    
    handles = guihandles(fh);
    
    val = get(handles.hstate,'Value');
    val = 3-val; %1->2 or 2->1
    set(handles.hstate,'Value',val);
    
    if val == 2   %if switching to snap state
        set(fh,'CurrentAxes',handles.ha);
        point = get(gca,'CurrentPoint');
        yfreeze = round(point(1,1));    %note the inverse convention for mouse and array (xy -> yx)
        xfreeze = round(point(1,2));
        xfreeze = max(1,min(xfreeze,nx));
        yfreeze = max(1,min(yfreeze,ny));
        
    end
    setappdata(fh,'xfreeze',xfreeze);
    setappdata(fh,'yfreeze',yfreeze);
    
    clear('point','val','point','xfreeze','yfreeze');
    clear('nx','ny','handles');
end

function change_Iminmax(hObj,eventdata,stack1,fh)
    
    handles = guihandles(fh);
    
    %replotting the main window
    plot_current_z_stack_single_channel(fh,stack1,1);
    
    %replotting the zoom window
    zoomwidth = str2num(get(handles.hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    plot_zoom_single_channel(x,y,1,zoomwidth,stack1,fh);
    
    clear('x','y','z','zoomwidth');
    clear('handles');
end

function global_auto_adjust(hglob_auto,eventdata,stack1,fh)
    
    handles = guihandles(fh);
    
    %updating the Imin/Imax info
    Imin = min(min(stack1)); 
    Imax = max(max(stack1));
    set(handles.hImin,'String',num2str(Imin)); 
    set(handles.hImax,'String',num2str(Imax));
    
    %replotting the main window
    plot_current_z_stack_single_channel(fh,stack1,1);
    
    %replotting the zoom window
    zoomwidth = str2num(get(handles.hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    
    plot_zoom_single_channel(x,y,1,zoomwidth,stack1,fh);
    
    clear('x','y','z','zoomwidth','Imin','Imax');
    clear('handles');
end

function viewlocaldata(src,eventdata,stack1,fh)
    
    nx = getappdata(fh,'nx');
    ny = getappdata(fh,'ny');
    
    handles = guihandles(fh);
    val = get(handles.hstate,'Value');          %this is the snap/grab mode
    
    %I replot things only if the grab mode is active
    if val == 1 

        %I get the current (x,y,z) 
        [x,y] = get_current_pointer_position(fh); 
        
        %I update the current display of the pointer info
        set(handles.hxpos,'String',num2str(x));
        set(handles.hypos,'String',num2str(y));

        if(x>=1 && y >=1 && x <= nx && y <= ny)
            Int = stack1(x,y);
            set(handles.hIval,'String',num2str(Int));

            %plotting the profiles and zoom window
            plot_local_data_single_channel(x,y,1,fh,stack1);
        end
    end
    clear('x','y','z','val','Int','nx','ny','handles');
end

function plot_local_data_single_channel(x,y,z,fh,stack1)
    
    handles = guihandles(fh);
    
    xrange = str2num(get(handles.hxplot_range,'String'));
    yrange = str2num(get(handles.hyplot_range,'String'));
    zoomwidth = str2num(get(handles.hzoom_width,'String'));

    %plotting the profiles and zoom window
    plot_xprofile(x,y,1,xrange,stack1,fh);
    plot_yprofile(x,y,1,yrange,stack1,fh);
    plot_zoom_single_channel(x,y,1,zoomwidth,stack1,fh);
    
    clear('handles','xrange','yrange','zrange','zoomwidth');
end

function plot_local_data_and_fit_single_channel_bgcorr(x,y,z,fh,stack1,p,Gaussout)
    
    handles = guihandles(fh);
    xrange = str2double(get(handles.hxplot_range,'String'));
    yrange = str2double(get(handles.hyplot_range,'String'));
    zoomwidth = str2double(get(handles.hzoom_width,'String'));

    %plotting the profiles and zoom window
    plot_xprofile_and_fit_bgcorr(x,y,1,xrange,stack1,fh,p,Gaussout);
    plot_yprofile_and_fit_bgcorr(x,y,1,yrange,stack1,fh,p,Gaussout);
    plot_zoom_single_channel(x,y,1,zoomwidth,stack1,fh);
    
    clear('xrange','yrange','zrange','zoomwidth','handles');
    
end

function plot_xprofile_and_fit_bgcorr(xarr,yarr,z,xrange,stack1,fh,p,Gaussout)

    handles = guihandles(fh);
    xh = handles.xh;
    nx = size(stack1,1);
    xprofile(1:nx) = stack1(1:nx,yarr);
    switch p.type
        case 'gaussian'
            xfit = Gaussout(2)+ Gaussout(1)*intensity_gaussian2D((1:nx)',yarr,...
                Gaussout(4),Gaussout(5),Gaussout(3),1,1);
        case 'integrated gaussian'
            I0 = intensity_integrated_gaussian2D(...
                1,1,0.5,0.5,Gaussout(3),1,1);
            xfit = Gaussout(2)+ Gaussout(1)*intensity_integrated_gaussian2D((1:nx)',yarr,...
                Gaussout(4),Gaussout(5),Gaussout(3),1,1)/I0;
    end
    xfit = xfit + Gaussout(8)*(1:nx)' + Gaussout(9)*yarr + Gaussout(10);
    
    if 2*xrange+1>nx
            set(fh,'CurrentAxes',handles.xh);
            Imin = min(xprofile(1:nx));
            Imax = max(xprofile(1:nx));
            plot(1:nx,xprofile(1:nx),[xarr xarr],[Imin Imax]); hold;
            plot(1:nx,xfit(1:nx),'r',[xarr xarr],[Imin Imax]); hold;        
    elseif (xarr - xrange >= 1) && (xarr + xrange <= nx)
            set(fh,'CurrentAxes',handles.xh);
            Imin = min(xprofile(xarr-xrange:xarr+xrange));
            Imax = max(xprofile(xarr-xrange:xarr+xrange));
            plot(xarr-xrange:xarr+xrange,xprofile(xarr-xrange:xarr+xrange),...
                [xarr xarr],[Imin Imax]); hold;
            plot(xarr-xrange:xarr+xrange,xfit(xarr-xrange:xarr+xrange),'r',...
                [xarr xarr],[Imin Imax]); hold;
    elseif xarr - xrange < 1 && xarr + xrange <= nx 
            set(fh,'CurrentAxes',handles.xh);
            x = xrange+1;
            Imin = min(xprofile(x-xrange:x+xrange));
            Imax = max(xprofile(x-xrange:x+xrange));
            plot(x-xrange:x+xrange,xprofile(x-xrange:x+xrange),...
                [xarr xarr],[Imin Imax]); hold;
            plot(x-xrange:x+xrange,xfit(x-xrange:x+xrange),'r',...
                [xarr xarr],[Imin Imax]); hold;        
    elseif xarr + xrange > nx && xarr - xrange >= 1
            set(fh,'CurrentAxes',handles.xh);
            x = nx-xrange;
            Imin = min(xprofile(x-xrange:x+xrange));
            Imax = max(xprofile(x-xrange:x+xrange));
            plot(x-xrange:x+xrange,xprofile(x-xrange:x+xrange),...
                [xarr xarr],[Imin Imax]); hold; 
            plot(x-xrange:x+xrange,xfit(x-xrange:x+xrange),'r',...
                [xarr xarr],[Imin Imax]); hold; 
    end
    set(xh,'Tag','xh');
    clear('nx','xprofile','Imin','Imax','x','pos','xfit');
    clear('handles');
end

function plot_yprofile_and_fit_bgcorr(xarr,yarr,z,yrange,stack1,fh,p,Gaussout)
    
    handles = guihandles(fh);
    yh = handles.yh;
    ny = size(stack1,2);
    yprofile(1:ny) = stack1(xarr,1:ny);
    
    switch p.type
        case 'gaussian'
            yfit = Gaussout(2)+ Gaussout(1)*intensity_gaussian2D(xarr,(1:ny)',...
                Gaussout(4),Gaussout(5),Gaussout(3),1,1);
        case 'integrated gaussian'
            I0 = intensity_integrated_gaussian2D(...
                1,1,0.5,0.5,Gaussout(3),1,1);
            yfit = Gaussout(2)+ Gaussout(1)*intensity_integrated_gaussian2D(xarr,(1:ny)',...
                Gaussout(4),Gaussout(5),Gaussout(3),1,1)/I0;
    end
    yfit = yfit + Gaussout(8)*xarr + Gaussout(9)*(1:ny)' + Gaussout(10);
    
    if 2*yrange+1>ny
            set(fh,'CurrentAxes',handles.yh);
            Imin = min(yprofile(1:ny));
            Imax = max(yprofile(1:ny));
            plot(1:ny,yprofile(1:ny),[yarr yarr],[Imin Imax]);hold;
            plot(1:ny,yfit(1:ny),'r',[yarr yarr],[Imin Imax]);hold;     
    elseif yarr - yrange >= 1 && yarr + yrange <= ny
            set(fh,'CurrentAxes',handles.yh);
            Imin = min(yprofile(yarr-yrange:yarr+yrange));
            Imax = max(yprofile(yarr-yrange:yarr+yrange));
            plot(yarr-yrange:yarr+yrange,yprofile(yarr-yrange:yarr+yrange),...
                [yarr yarr],[Imin Imax]);hold;
            plot(yarr-yrange:yarr+yrange,yfit(yarr-yrange:yarr+yrange),'r',...
                [yarr yarr],[Imin Imax]);hold;
    elseif yarr - yrange < 1  && yarr + yrange <= ny
            set(fh,'CurrentAxes',handles.yh);
            y = yrange+1;
            Imin = min(yprofile(y-yrange:y+yrange));
            Imax = max(yprofile(y-yrange:y+yrange));
            plot(y-yrange:y+yrange,yprofile(y-yrange:y+yrange),...
                [yarr yarr],[Imin Imax]);hold;
            plot(y-yrange:y+yrange,yfit(y-yrange:y+yrange),'r',...
                [yarr yarr],[Imin Imax]);hold;
    elseif yarr + yrange > ny && yarr - yrange >= 1
            set(fh,'CurrentAxes',handles.yh);
            y = ny-yrange;
            Imin = min(yprofile(y-yrange:y+yrange));
            Imax = max(yprofile(y-yrange:y+yrange));
            plot(y-yrange:y+yrange,yprofile(y-yrange:y+yrange),...
                [yarr yarr],[Imin Imax]); hold;
            plot(y-yrange:y+yrange,yfit(y-yrange:y+yrange),'r',...
                [yarr yarr],[Imin Imax]); hold;
    end 
    set(yh,'Tag','yh');
    clear('ny','yprofile','Imin','Imax','y','pos','yfit');
    clear('xarr','yarr','handles');
end

function [xarr yarr] = get_current_pointer_position(fh)

    xfreeze = getappdata(fh,'xfreeze');
    yfreeze = getappdata(fh,'yfreeze');
    nx = getappdata(fh,'nx');
    ny = getappdata(fh,'ny');
    
    handles = guihandles(fh);
    
    val = get(handles.hstate,'Value');%this is the snap/grab mode for the local data panel

    if val == 1         %if grab mode is activated, I replot the local data at some arbitrary location
            set(fh,'CurrentAxes',handles.ha);
            point = get(gca,'CurrentPoint');
            xpic = round(point(1,1));   
            ypic = round(point(1,2));
            xarr = ypic;    
            yarr = xpic;    %pic and array have different xy conventions
            xarr = max(xarr,1);  
            xarr = min(xarr,nx);
            yarr = max(yarr,1);  
            yarr = min(yarr,ny);
    else
            xarr = xfreeze; 
            yarr = yfreeze;
    end

    clear('xfreeze','yfreeze','nx','ny','handles','xpic','ypic','val');
end

function plot_xprofile(xarr,yarr,z,xrange,stack1,fh)
    
    handles = guihandles(fh);
    nx = size(stack1,1);
    xprofile(1:nx) = stack1(1:nx,yarr);
    xh = handles.xh;
    if 2*xrange+1>nx
            set(fh,'CurrentAxes',handles.xh);
            Imin = min(xprofile(1:nx));
            Imax = max(xprofile(1:nx));
            plot(1:nx,xprofile(1:nx),[xarr xarr],[Imin Imax]); 
    elseif (xarr - xrange >= 1) && (xarr + xrange <= nx)
            set(fh,'CurrentAxes',handles.xh);
            Imin = min(xprofile(xarr-xrange:xarr+xrange));
            Imax = max(xprofile(xarr-xrange:xarr+xrange));
            plot(xarr-xrange:xarr+xrange,xprofile(xarr-xrange:xarr+xrange),...
                [xarr xarr],[Imin Imax]);
    elseif xarr - xrange < 1 && xarr + xrange <= nx 
            set(fh,'CurrentAxes',handles.xh);
            x = xrange+1;
            Imin = min(xprofile(x-xrange:x+xrange));
            Imax = max(xprofile(x-xrange:x+xrange));
            plot(x-xrange:x+xrange,xprofile(x-xrange:x+xrange),...
                [xarr xarr],[Imin Imax]);
    elseif xarr + xrange > nx && xarr - xrange >= 1
            set(fh,'CurrentAxes',handles.xh);
            x = nx-xrange;
            Imin = min(xprofile(x-xrange:x+xrange));
            Imax = max(xprofile(x-xrange:x+xrange));
            plot(x-xrange:x+xrange,xprofile(x-xrange:x+xrange),...
                [xarr xarr],[Imin Imax]);   
    end 
    set(xh,'Tag','xh');
    clear('nx','xprofile','Imin','Imax','x','handles','xh');
end

function plot_yprofile(xarr,yarr,z,yrange,stack1,fh)

    handles = guihandles(fh);
    yh = handles.yh;
    ny = size(stack1,2);
    yprofile(1:ny) = stack1(xarr,1:ny);

    if 2*yrange+1>ny
            set(fh,'CurrentAxes',handles.yh);
            Imin = min(yprofile(1:ny));
            Imax = max(yprofile(1:ny));
            plot(1:ny,yprofile(1:ny),[yarr yarr],[Imin Imax]);
    elseif yarr - yrange >= 1 && yarr + yrange <= ny
            set(fh,'CurrentAxes',handles.yh);
            Imin = min(yprofile(yarr-yrange:yarr+yrange));
            Imax = max(yprofile(yarr-yrange:yarr+yrange));
            plot(yarr-yrange:yarr+yrange,yprofile(yarr-yrange:yarr+yrange),...
                [yarr yarr],[Imin Imax]);
    elseif yarr - yrange < 1  && yarr + yrange <= ny
            set(fh,'CurrentAxes',handles.yh);
            y = yrange+1;
            Imin = min(yprofile(y-yrange:y+yrange));
            Imax = max(yprofile(y-yrange:y+yrange));
            plot(y-yrange:y+yrange,yprofile(y-yrange:y+yrange),...
                [yarr yarr],[Imin Imax]);
    elseif yarr + yrange > ny && yarr - yrange >= 1
            set(fh,'CurrentAxes',handles.yh);
            y = ny-yrange;
            Imin = min(yprofile(y-yrange:y+yrange));
            Imax = max(yprofile(y-yrange:y+yrange));
            plot(y-yrange:y+yrange,yprofile(y-yrange:y+yrange),...
                [yarr yarr],[Imin Imax]);
    end 
    set(yh,'Tag','yh');
    clear('ny','yprofile','Imin','Imax','y','handles','yh');
end

function I = gauss_integrated2D_position_list_bgcorr(Coeffs, pos)
%everything is in pixel units here
%origin of positions at the corner on the image 
%(i.e. centers of pixels have half integer vox_size values)
%pos is an array of positions [x y z]

Itot = Coeffs(6); 
bg = Coeffs(2);
sxy=Coeffs(3); 
xc=Coeffs(4); 
yc=Coeffs(5);

r = pos(:,1);
c = pos(:,2);

diffx1 =  (double(r-1)) - xc;
diffx1 = diffx1 ./ ( sqrt(2) * sxy );
diffx2 =  double(r) - xc;
diffx2 = diffx2 ./ ( sqrt(2) * sxy );

intensity1 = abs( erf( diffx1) - erf(diffx2) );

diffy1 =  (double(c-1)) - yc;
diffy1 = diffy1 ./ ( sqrt(2) * sxy );
diffy2 =  double(c) - yc;
diffy2 = diffy2 ./ ( sqrt(2) * sxy );

intensity2 = abs( erf( diffy1) - erf(diffy2) );

intensity = Itot.*intensity1.*intensity2/4.0;
I = intensity + bg;

%I = I + Coeffs(8)*(double(r) - Coeffs(11)) + Coeffs(9)*(double(c) - Coeffs(12)) + Coeffs(10);
I = I + Coeffs(8)*double(r)  + Coeffs(9)*double(c) + Coeffs(10);

clear('Imax','r','c','h','diffx1','diffx2','diffy1','diffy2','diffz1','diffz2',...
    'intensity2','intensity3','intensity','bg','sxy','sz','xc','yc','zc');
end

function fh = set_figure_and_panel_for_explore_img(stack1,stack1name,fh)


%%%%%%%%%%%%%%%%%%%%% variables initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

strXpos = '1';     % default x value of the current pixel
strYpos = '1';     % default y value of the current pixel
strIval = num2str(stack1(1,1,1));   % value of the intensity @ the current pixel

%values of global min max and median
Imed = median_im(double(stack1));
Isd = std_im(double(stack1));
Imin = min(min(stack1));
Imax = max(max(stack1));

strImin = num2str(Imin,'%10.4f');
strImax = num2str(Imax,'%10.4f');
strImed = num2str(Imed,'%10.4f');
strIsd = num2str(Isd,'%10.4f');

%default values of the ranges of the 3 axial views
xpix = num2str(50);
ypix = num2str(50);

%default value of the size of the zoom window.
startzoomwidth = num2str(50);

%% %%%%%%%%%%%%%%%%%%%%%%% main figure, plots and title %%%%%%%%%%%%%%%%%%%%%%%
 

%set(fh,'Toolbar','figure');
set(fh,'Name',stack1name);              
axes('Parent',fh,'Units','characters','Position',[10,6.15,102.4,39.4],'Tag','ha');     %main figure
        
axes('Parent',fh,'Units','characters','Position',[120 45.4 66 7.7],'Tag','xh'); 
axes('Parent',fh,'Units','characters','Position',[120 33.8 66 7.7],'Tag','yh'); 

%% %%%%%%%%%%%%%%%%%%%%%%%%% settings for the axial views %%%%%%%%%%%%%%%%
uicontrol('Units','characters',...
                'Style','text','String','profile along x',...
                'Position',[144,53.3,20,1.2]);
uicontrol('Units','characters',...
                'Style','text','String','profile along y',...
                'Position',[144,41.8,20,1.2]);

uicontrol('Units','characters',...
                'Style','text','String','number of voxels to display',...
                'Position',[188,49.2,10,3.6]); 
uicontrol('Units','characters',...
                'Style','text','String','number of voxels to display',...
                'Position',[188,37.7,10,3.6]); 

uicontrol('Units','characters',...
                'Tag','hxplot_range',...
                'Style','edit','String',xpix,...
                'Position',[190,47.7,6,1.2]);                  
uicontrol('Units','characters',...
                'Tag','hyplot_range',...
                'Style','edit','String',ypix,...
                'Position',[190,36.2,6,1.2]);            


%% %%%%%%%%%%%%%%%%%%%%%%%%% zoom plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
axes('Parent',fh,'Units','characters','Position',[144 1.9 46 17.7],'Tag','zoomh');             

uicontrol('Units','characters',...
                'Style','text','String','window size',...
                'Position',[192,12,8,2.3]);
            
uicontrol('Units','characters',...
                'Tag','hzoom_width',...
                'Style','edit','String',startzoomwidth,...
                'Position',[192,10.4,6,1.2]); 

%% %%%%%%%%%%%%%%%%%%%%%%%% fit gaussian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uicontrol('Units','characters',...
                'Tag','hgauss_fit',...
                'Style','pushbutton','String','local gaussian fit',...
                'Position',[116,16.7,22,2]);

uicontrol('Units','characters',...
                'Style','text','String','x max',...
                'Position',[116,15.4,8,1.2]);
uicontrol('Units','characters',...
                'Tag','hx0',...
                'Style','text','String','undefined',...
                'Position',[126,15.4,12,1.2]);
uicontrol('Units','characters',...
                'Style','text','String','y max',...
                'Position',[116,14.2,8,1.2]);
uicontrol('Units','characters',...
                'Tag','hy0',...
                'Style','text','String','undefined',...
                'Position',[126,14.2,12,1.2]);            
uicontrol('Units','characters',...
                'Style','text','String','Int. max',...
                'HorizontalAlignment','left',...
                'Position',[116,11.5,14,1.2]);
uicontrol('Units','characters',...
                'Tag','hI0',...
                'Style','text','String','undefined',...
                'Position',[126,11.5,12,1.2]); 
uicontrol('Units','characters',...
                'Style','text','String','background',...
                'HorizontalAlignment','left',...
                'Position',[116,10.4,16,1.2]);
uicontrol('Units','characters',...
                'Tag','hbg0',...
                'Style','text','String','undefined',...
                'Position',[132,10.4,6,1.2]); 
uicontrol('Units','characters',...
                'Style','text','String','s_xy',...
                'Position',[116,8.8,8,1.2]);
uicontrol('Units','characters',...
                 'Tag','hsxy0',...
                'Style','text','String','undefined',...
                'Position',[126,8.8,12,1.2]);              
            
uicontrol('Units','characters',...
                'Tag','hrecord',...
                'Style','PushButton','String','Record Fit Results',...
                'Position',[116,5.5,22,2]);                       
uicontrol('Units','characters',...
                'Style','text','String','Fits Recorded',...
                'HorizontalAlignment','left',...
                'Position',[116,4.1,14,1.2]);
uicontrol('Units','characters',...
                'Tag','hFitsRecorded',...
                'Style','text','String','0',...
                'Position',[132,4.1,6,1.2]);            
                 
%% %%%%%%%%%%%%%%%%%%%%%%%%% contrast panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ph = uipanel(fh,'Title','Contrast','Units','characters',...
             'Position',[30 46.2 83 8.5],'TitlePosition','centertop');


            
%info on whole stack                        
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','Image',...
                'Position',[10,2.3,18,1.2]);        

uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','min',...
                'Position',[0,1.2,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','median',...
                'Position',[13,1.2,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','max',...
                'Position',[24,1.2,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','SD',...
                'Position',[39,1.2,8,1.2]);
            
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hglobal_Imin_Info_val',...
                'Style','text','String',strImin,...
                'Position',[0,1,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hglobal_Imed_Info_val',...
                'Style','text','String',strImed,...
                'Position',[13,1,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hglobal_Imax_Info_val',...
                'Style','text','String',strImax,...
                'Position',[24,1,8,1.2]);                       
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hglobal_Isd_Info_val',...
                'Style','text','String',strIsd,...
                'Position',[39,1,8,1.2]);             
            
%contrast settings           


uicontrol('Parent',ph,'Units','characters',...
                'Tag','hglob_auto',...
                'Style','pushbutton',...
                'String','auto adjust',...
                'Units','characters','Position',[52,5,30,1.9]); 
            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','min Int.',...
                'Position',[52,1.3,12,1.2]);
       
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hImin',...
                'Style','edit','String',strImin,...
                'Position',[52,0.4,12,1.2]);         

uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','max Int.',...
                'Position',[66,1.3,12,1.2]);

uicontrol('Parent',ph,'Units','characters',...
                'Tag','hImax',...
                'Style','edit','String',strImax,...
                'Position',[66,0.4,12,1.2]);                
       
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% pointer info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uicontrol('Style','text','String','x',...
            'Units','characters',...
           'Position',[10,48.5,6,1.2]);

uicontrol('style','list','String',...
{'<HTML><FONT COLOR=00FF00>grab</FONT><HTML>',...
'<HTML><FONT COLOR=FF0000>snap</FONT><HTML>'},...
'Units','characters',...
'Tag','hstate',...
'Position',[120,53.3,12,1.2]);       

uicontrol('Style','text','String',strXpos,...
            'Units','characters',...
            'Tag','hxpos',...
           'Position',[14,48.5,12,1.2]);       
       
uicontrol('Style','text','String','y',...
           'Units','characters',...
           'Position',[10,47.3,6,1.2]);

uicontrol('Style','text','String',strYpos,...
            'Tag','hypos',...
            'Units','characters',...
           'Position',[14,47.3,12,1.2]);       
       
uicontrol('Style','text','String','Int',...
            'Units','characters',...
           'Position',[10,46.2,6,1.2]);

uicontrol('Style','text','String',strIval,...
            'Units','characters',...
            'Tag','hIval',...
           'Position',[14,46.2,12,1.2]);


uicontrol(fh,'Style','pushbutton',...
                    'Tag','hclose',...
                    'Units','characters',...
                    'String','done',...
                    'Position',[116 0.4 22 3]);
                             
clear('nz','strXpos','strYpos','strZpos',...
    'Imed','Isd','Imin','Imax',...
    'strImin','strImax','strImed','strIsd',...
    'curImin','curImax','curImed','curIsd',...
    'strcurImin','strcurImax','strcurImed','strcurIsd',...
    'xpix','ypix','zpix',...
    'startzoomwidth','ph');

clear('stack1','stack1name');
end
