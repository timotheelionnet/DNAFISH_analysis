function params = explore_stack_2channels2(params)

%% creating figure
fh = figure(...
              'Units','characters',...
              'MenuBar','none',...
              'Toolbar','none',...
              'NumberTitle','off',...
              'Position',[20 3 200 55],...
              'Visible','off'); 
          
%% parsing arguments 
%checking that we have at least one 3D stack or load one. this is the data (green channel)
 
stack1name = inputname(1);
stack1 = params.data;
params = rmfield(params,'data');
stack2 = params.smooth;
params = rmfield(params,'smooth');

if ndims(stack1) ~= 3 || size(stack1,3)==1
    fprintf(['your stack size is', num2str(size(stack1)),' it should be 3D.\n']);
    return;
end

[nx ny nz] = size(stack1);
setappdata(fh,'nx',nx);
setappdata(fh,'ny',ny);
setappdata(fh,'nz',nz);
xfreeze = round(nx/2); 
yfreeze = round(ny/2); 
setappdata(fh,'xfreeze',xfreeze);
setappdata(fh,'yfreeze',yfreeze);
params.thresh.sd = std_stack(stack2);
setappdata(fh,'thresh',params.thresh);


%storing each z plane separately so that I can save memory when displaying
%them
nz = size(stack2,3);
for i=1:nz
    eval(['st',num2str(i),' = stack2(:,:,i);']);
    eval(['setappdata(fh,''st',num2str(i),''',st',num2str(i),');']);
    eval(['clear(''st',num2str(i),''');']);
end
clear('stack2');

fh = set_figure_and_panel_for_explore_stack_set_threshold(stack1,stack1name,fh,params);
handles = guihandles(fh);


change_overlay_channel_clean2(fh);


%%%%%%%%%%%%%%%%%%%%%%%%%%% set all the callback functions %%%%%%%%%%%%%%%%%            
set(fh,'WindowButtonMotionFcn', {@viewlocaldata,stack1,fh});
set(fh,'WindowButtonUpFcn', {@change_state,fh});   
set(handles.hclose,'Callback', {@clean_close,fh});   
set(handles.zslider,'Callback',{@slider_change,stack1,fh});
set(handles.hzpos,'Callback',{@zpos_change,stack1,fh});            
set(handles.hImin,'Callback',{@change_Iminmax,stack1,fh});
set(handles.hImax,'Callback',{@change_Iminmax,stack1,fh});
set(handles.hslice_auto,'Callback',{@slice_auto_adjust,stack1,fh});
set(handles.hglob_auto,'Callback',{@global_auto_adjust,stack1,fh});
set(handles.hxplot_range,'Callback',{@change_local_range,stack1,fh});
set(handles.hyplot_range,'Callback',{@change_local_range,stack1,fh});
set(handles.hzplot_range,'Callback',{@change_local_range,stack1,fh});
set(handles.hoverlay,'Callback',{@toggle_overlay,fh,stack1});
set(handles.hthresh_level,'Callback',{@change_thresh,fh,stack1});
set(handles.hthresh_units,'SelectionChangeFcn',{@switch_thresh_units,fh,stack1});

set(fh,'SelectionType','alt');
set(fh,'CurrentAxes',handles.ha);
    
plot_current_z_stack_two_channels(fh,stack1,1);
set(fh,'Visible','on');

uiwait;

params.thresh = getappdata(fh,'thresh');
params.data = stack1;
clear('stack1');

%mem00 = memory
%params.smooth = zeros(size(params.data));
%mem0 = memory
for i=1:nz
    %disp(num2str(i))
    %mem = memory
    %eval(['params.smooth(:,:,i) = getappdata(fh,''st',num2str(i),''');']);
    %mem = memory
    %eval(['rmappdata(fh,''st',num2str(i),''');']);
end

close(fh);
drawnow;
%% housekeeping
clear('nx','ny','nz','stack1','stack1name','xfreeze','yfreeze');
clear('fh','ha','handles');

end

%%%%%%%%%%%%%%%%%%%%%%%%%% Callback functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function change_thresh(src,eventdata,fh,stack1)
    handles = guihandles(fh);
    
    %updating the threshold value
    thresh = getappdata(fh,'thresh');
    thresh.level = str2double(get(handles.hthresh_level,'String'));
    setappdata(fh,'thresh',thresh);
    
    change_overlay_channel_clean2(fh);
    npix_txt = num2str(getappdata(fh,'npix_above_thresh'));
    set(handles.hnpix_above_threshold,'String',npix_txt);
    
    %replotting the local data (zoom + axial views)
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    z = str2double(get(handles.hzpos,'String'));
    z = round(z);
    
    plot_local_data_two_channels(x,y,z,fh,stack1);
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,stack1,z);
    
    %cleaning the global link 
    clear('thresh','stack2','xfreeze','yfreeze','handles','thresh',...
        'x','y','z','npix_txt','npix_above_threshold');
end

function switch_thresh_units(src,eventdata,fh,stack1)
    handles = guihandles(fh);
    
    %updating the threshold value
    thresh = getappdata(fh,'thresh');
    
    str_units = get(get(handles.hthresh_units,'SelectedObject'),'Tag');
    if strcmp(str_units,'absolute') 
        thresh.units = 'absolute';
    elseif strcmp(str_units,'SD') 
        thresh.units = 'SD';
    elseif strcmp(str_units,'automatic') 
        thresh.units = 'automatic';
    end
    
    setappdata(fh,'thresh',thresh);
    
    change_overlay_channel_clean2(fh);
    npix_txt = num2str(getappdata(fh,'npix_above_thresh'));
    set(handles.hnpix_above_threshold,'String',npix_txt);
    
    %replotting the local data (zoom + axial views)
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    z = str2double(get(handles.hzpos,'String'));
    z = round(z);
    
    plot_local_data_two_channels(x,y,z,fh,stack1);
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,stack1,z);
    
    %cleaning the global link 
    clear('thresh','xfreeze','yfreeze','handles','thresh',...
        'stdev','threshInt','npix_above_threshold','npix_txt','x','y','z');
end

function change_overlay_channel_clean2(fh)
    thresh = getappdata(fh,'thresh');
    
    if strcmp(thresh.units,'absolute')
        threshInt = thresh.level;
    elseif strcmp(thresh.units,'SD')
        threshInt = thresh.level*thresh.sd;
    elseif strcmp(thresh.units,'automatic')
        threshInt = thresh.thr_min + thresh.level*(thresh.thr_opt - thresh.thr_min);
    end
    
    tmp = 0;
    for i=1:getappdata(fh,'nz')
        eval(['st = getappdata(fh,''st',num2str(i),''');']);
        tmp = tmp + sum(sum(st > threshInt));
    end
    
    setappdata(fh,'npix_above_thresh',tmp );
    
    clear('thresh','threshInt','tmp');
end

function toggle_overlay(src,eventdata,fh,stack1)
    handles = guihandles(fh);   
    [~, ny] = size(stack1);
    
    z = str2double(get(handles.hzpos,'String'));
    z = round(z);
    
    %switching the state of the button
    is_overlay_on = 1 - get(handles.hoverlay,'Value');
    if is_overlay_on == 1 
        set(handles.hoverlay,'String','hide overlay');
    else 
        set(handles.hoverlay,'String','show overlay');
    end
    
    %replotting the zoom window
    zoomwidth = str2double(get(handles.hzoom_width,'String'));
    
    val = get(handles.hstate,'Value');  %this is the snap/grab mode for the local data panel
    if val == 1                 %if grab mode is activated, I replot the zoom data at some arbitrary location
        plot_zoom_two_channels(1,ny,z,zoomwidth,stack1,fh);
          
    else                %if snap mode is activated, I plot the local data at the same location
        plot_zoom_two_channels(getappdata(fh,'xfreeze'),getappdata(fh,'yfreeze'),z,zoomwidth,stack1,fh);
    end
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,stack1,z);
    
    clear('xfreeze','yfreeze','stack2','handles',...
        'nx','ny','nz','is_overlay_on','z','val','zoomwidth');
end

function plot_current_z_stack_two_channels(fh,stack1,z)

thresh = getappdata(fh,'thresh');
if strcmp(thresh.units,'absolute')
    threshInt = thresh.level;
elseif strcmp(thresh.units,'SD')
    threshInt = thresh.level*thresh.sd;
elseif strcmp(thresh.units,'automatic')
    threshInt = thresh.thr_min + thresh.level*(thresh.thr_opt - thresh.thr_min);
end

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
    
    %getting the plane of interest
    eval(['stack2 = getappdata(fh,''st',num2str(z),''') - threshInt > 0;']);

%% filling data 
    cur_slice = zeros(nx,ny,3);
    
    %red channel: overlay or nothing
    is_overlay_on = 1 - get(handles.hoverlay,'Value');
    if is_overlay_on
        cur_slice(1:nx,1:ny,1) = stack2(1:nx,1:ny);
    else
        cur_slice(1:nx,1:ny,1) = 0;
    end
    
    %green channel: data
    if(Imin >= Imax)
        cur_slice(1:nx,1:ny,2) = 0;
    else
        cur_slice(1:nx,1:ny,2) = min( max((double(stack1(1:nx,1:ny))-Imin)/(Imax-Imin),0) , 1);
    end
    
    %nothing in the blue channel
    cur_slice(1:nx,1:ny,3) = 0; 
    
%%    
    ha = handles.ha;
    set(fh,'CurrentAxes',handles.ha);
    
    %ha = imagesc(cur_slice,'XLim',[1,max(nx,ny)],'YLim',[1,max(nx,ny)]); colormap(gray);
    imagesc(cur_slice); 
    xlim([1,max(nx,ny)]);
    ylim([1,max(nx,ny)]);
    set(ha,'Tag','ha');
    
    %set(ha,'XLim',[1,max(nx,ny)] );
    %set(ha,'YLim',[1,max(nx,ny)] );

clear('cur_slice','nx','ny','nz','Imin','Imax','newXLim','newYLim','stack1');
clear('stack2','handles','ha','thresh','threshInt','is_overlay_on');
end

function plot_zoom_two_channels(x,y,z,width,stack1,fh)
 %inputs a 3D or 2D image and plots a locally  zoomed version  
    handles = guihandles(fh);
    set(fh,'CurrentAxes',handles.zoomh);
    zoomh = handles.zoomh;
    thresh = getappdata(fh,'thresh');
    if strcmp(thresh.units,'absolute')
        threshInt = thresh.level;
    elseif strcmp(thresh.units,'SD')
        threshInt = thresh.level*thresh.sd;
    elseif strcmp(thresh.units,'automatic')
        threshInt = thresh.thr_min + thresh.level*(thresh.thr_opt - thresh.thr_min);
    end
    

    
    if ndims(stack1) == 3
        [nx ny nz] = size(stack1);
        if z<=0 || z >nz
            disp('z out of bounds in plot_current_z_stack!');
        end
        stack1 = stack1(:,:,z); 
        
    elseif ndims(stack1) ==2
        [nx ny] = size(stack1);
    end
    
    %getting the plane of interest
    eval(['stack2 = getappdata(fh,''st',num2str(z),''') - threshInt > 0;']);
    
    Imin = str2num(get(handles.hImin,'String'));
    Imax = str2num(get(handles.hImax,'String'));

%% computing the value of the local slice 
    
    cur_slice = zeros(nx,ny,3);

    %red channel: overlay or nothing
    is_overlay_on = 1 - get(handles.hoverlay,'Value');
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
    
%% chosing the window corresponding to the cursor position
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
    set(zoomh,'Tag','zoomh');
    clear('cur_slice','nx','ny','nz','Imin','Imax','is_overlay_on','stack1','stack2');
    clear('xx','yy','handles','thresh','threshInt','zoomh');
    clear('x','y','z','width','zoomh','stack1','stack2','fh','hImin','hImax','hoverlay');
end

function change_local_range(src,eventdata,stack1,fh)
    handles = guihandles(fh);
    z = get(handles.zslider,'Value');
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    plot_local_data_two_channels(x,y,z,fh,stack1);
    clear('x','y','z','handles');
end

function clean_close(src,eventdata,fh)
    uiresume;
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

function slider_change(zslider,eventdata,stack1,fh)
    
    xfreeze = getappdata(fh,'xfreeze');
    yfreeze = getappdata(fh,'yfreeze');
    nx = getappdata(fh,'nx');
    ny = getappdata(fh,'ny');
    
    handles = guihandles(fh);
    
    %updating the z-slider and the z info window
    z = get(handles.zslider,'Value');
    z = round(z);
    set(handles.zslider,'Value',z);
    set(handles.hzpos,'String',num2str(z));
    
    %updating the Imin/Imax info
    temp = reshape(stack1(1:nx,1:ny,z),1,nx*ny);
    curImin = min(temp); curImax = max(temp); curImed = median(double(temp)); curIsd = std(double(temp));
    clear('temp');
    set(handles.hslice_Imin_Info_val,'String',num2str(curImin));
    set(handles.hslice_Imed_Info_val,'String',num2str(curImed));
    set(handles.hslice_Imax_Info_val,'String',num2str(curImax)); 
    set(handles.hslice_Isd_Info_val,'String',num2str(curIsd)); 
    
    %replotting the local data (zoom + axial views)
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    
    Int = stack1(x,y,z);
    set(handles.hIval,'String',num2str(Int));

    plot_local_data_two_channels(x,y,z,fh,stack1);
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,stack1,z);  
    
    clear('x','y','z','Int','curImin','curImax','curImed','curIsd');
    clear('nx','ny','xfreeze','yfreeze','handles');
end

function zpos_change(hzpos,eventdata,stack1,fh)

    xfreeze = getappdata(fh,'xfreeze');
    yfreeze = getappdata(fh,'yfreeze');
    nx = getappdata(fh,'nx');
    ny = getappdata(fh,'ny');
    nz = getappdata(fh,'nz');
    
    handles = guihandles(fh);
    
    %updating the z-slider and the z info window
    z = str2num(get(handles.hzpos,'String'));
    z = round(z);
    z = max(min(z,nz) , 1);
    set(handles.hzpos,'String',num2str(z));
    set(handles.zslider,'Value',z);
    
    %updating the Imin/Imax info
    temp = reshape(stack1(1:nx,1:ny,2),1,nx*ny);
    curImin = min(temp); curImax = max(temp); curImed = median(double(temp)); curIsd = std(double(temp));
    clear('temp');
    set(handles.hslice_Imin_Info_val,'String',num2str(curImin));
    set(handles.hslice_Imed_Info_val,'String',num2str(curImed));
    set(handles.hslice_Imax_Info_val,'String',num2str(curImax));
    set(handles.hslice_Isd_Info_val,'String',num2str(curIsd)); 
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,stack1,z);
    
    %replotting the local data (zoom + axial views)
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    Int = stack1(x,y,z);
    set(handles.hIval,'String',num2str(Int));

    plot_local_data_two_channels(x,y,z,fh,stack1);
    
    clear('x','y','z','Int','curImin','curImax','curImed','curIsd'); 
    clear('nx','ny','xfreeze','yfreeze','handles');
end

function change_Iminmax(hObj,eventdata,stack1,fh)
    
    handles = guihandles(fh);
    z = get(handles.zslider,'Value');
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,stack1,z);
    
    %replotting the zoom window
    zoomwidth = str2num(get(handles.hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    plot_zoom_two_channels(x,y,z,zoomwidth,stack1,fh);
    
    clear('x','y','z','zoomwidth');
    clear('handles');
end

function slice_auto_adjust(hslice_auto,eventdata,stack1,fh)
    
    handles = guihandles(fh);
    z = get(handles.zslider,'Value');
    
    %updating Imin / Imax info
    Imin = min(min(stack1(:,:,z))); 
    Imax = max(max(stack1(:,:,z)));
    set(handles.hImin,'String',num2str(Imin)); 
    set(handles.hImax,'String',num2str(Imax));
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,stack1,z);
    
    %replotting the zoom window
    zoomwidth = str2num(get(handles.hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    
    plot_zoom_two_channels(x,y,z,zoomwidth,stack1,fh);
    clear('x','y','z','zoomwidth','Imin','Imax');
    clear('handles');
end

function global_auto_adjust(hglob_auto,eventdata,stack1,fh)
    
    handles = guihandles(fh);
    z = get(handles.zslider,'Value');
    
    %updating the Imin/Imax info
    Imin = min(min(min(stack1))); 
    Imax = max(max(max(stack1)));
    set(handles.hImin,'String',num2str(Imin)); 
    set(handles.hImax,'String',num2str(Imax));
    
    %replotting the main window
    plot_current_z_stack_two_channels(fh,stack1,z);
    
    %replotting the zoom window
    zoomwidth = str2num(get(handles.hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh); %I get the current (x,y) of the mouse
    
    plot_zoom_two_channels(x,y,z,zoomwidth,stack1,fh);
    
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
        z = get(handles.zslider,'Value');
        z = round(z);

        %I update the current display of the pointer info
        set(handles.hxpos,'String',num2str(x));
        set(handles.hypos,'String',num2str(y));

        if(x>=1 && y >=1 && x <= nx && y <= ny)
            Int = stack1(x,y,z);
            set(handles.hIval,'String',num2str(Int));

            %plotting the profiles and zoom window
            plot_local_data_two_channels(x,y,z,fh,stack1);
        end
    end
    clear('x','y','z','val','Int','nx','ny','handles');
end

function plot_local_data_two_channels(x,y,z,fh,stack1)
   
    handles = guihandles(fh);
    
    xrange = str2num(get(handles.hxplot_range,'String'));
    yrange = str2num(get(handles.hyplot_range,'String'));
    zrange = str2num(get(handles.hzplot_range,'String'));
    zoomwidth = str2num(get(handles.hzoom_width,'String'));

    %plotting the profiles and zoom window
    plot_xprofile(x,y,z,xrange,stack1,fh);
    plot_yprofile(x,y,z,yrange,stack1,fh);
    plot_zprofile(x,y,z,zrange,stack1,fh);
    plot_zoom_two_channels(x,y,z,zoomwidth,stack1,fh);
    
    clear('handles','xrange','yrange','zrange','zoomwidth');
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
    xprofile(1:nx) = stack1(1:nx,yarr,z);
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
    yprofile(1:ny) = stack1(xarr,1:ny,z);

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

function plot_zprofile(xarr,yarr,z,zrange,stack1,fh)
    
    handles = guihandles(fh);
    zh = handles.zh;
    
    nz = size(stack1,3);
    zprofile(1:nz) = stack1(xarr,yarr,1:nz);

    if 2*zrange+1 > nz
            set(fh,'CurrentAxes',handles.zh);
            Imin = min(zprofile(1:nz));
            Imax = max(zprofile(1:nz));
            plot(1:nz,zprofile(1:nz),[z z],[Imin Imax]);
    elseif z - zrange >= 1 && z + zrange <= nz
            set(fh,'CurrentAxes',handles.zh);
            Imin = min(zprofile(z-zrange:z+zrange));
            Imax = max(zprofile(z-zrange:z+zrange));
            plot(z-zrange:z+zrange,zprofile(z-zrange:z+zrange),[z z],[Imin Imax]);
    elseif z - zrange < 1 && z + zrange <= nz  
            set(fh,'CurrentAxes',handles.zh);
            zz = zrange+1;
            Imin = min(zprofile(zz-zrange:zz+zrange));
            Imax = max(zprofile(zz-zrange:zz+zrange));
            plot(zz-zrange:zz+zrange,zprofile(zz-zrange:zz+zrange),...
                [z z],[Imin Imax]);
    elseif z + zrange > nz && z - zrange >= 1
            set(fh,'CurrentAxes',handles.zh);
            zz = nz-zrange;
            Imin = min(zprofile(zz-zrange:zz+zrange));
            Imax = max(zprofile(zz-zrange:zz+zrange));
            plot(zz-zrange:zz+zrange,zprofile(zz-zrange:zz+zrange),...
                [z z],[Imin Imax]);
    end 
    set(zh,'Tag','zh');
    clear('nz','zprofile','Imin','Imax','z','zz','handles','zh');
end

function fh = set_figure_and_panel_for_explore_stack_set_threshold(stack1,stack1name,fh,params)


%%%%%%%%%%%%%%%%%%%%% variables initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,nz] = size(stack1);

strXpos = '1';     % default x value of the current pixel
strYpos = '1';     % default y value of the current pixel
strIval = num2str(stack1(1,1,1));   % value of the intensity @ the current pixel
strZpos = '1';     % the default value for the z stack displayed

%values of global min max and median
Imed = median(double(stack1(:)));
Isd = std(double(stack1(:)));
Imin = min(min(min(stack1)));
Imax = max(max(max(stack1)));

strImin = num2str(Imin,'%10.4f');
strImax = num2str(Imax,'%10.4f');
strImed = num2str(Imed,'%10.4f');
strIsd = num2str(Isd,'%10.4f');

%values of min max and median of current slice (z=1 when you open the viewer)
curImin = min(min(stack1(:,:,1)));
curImax = max(max(stack1(:,:,1)));
curImed = median_im(double(stack1(:,:,1)));
curIsd = std_im(stack1(:,:,1));

strcurImin = num2str(curImin,'%10.4f');
strcurImax = num2str(curImax,'%10.4f');
strcurImed = num2str(curImed,'%10.4f');
strcurIsd = num2str(curIsd,'%10.4f');

%default values of the ranges of the 3 axial views
xpix = num2str(50);
ypix = num2str(50);
zpix = num2str(50);

%default value of the size of the zoom window.
startzoomwidth = num2str(50);

%% %%%%%%%%%%%%%%%%%%%%%%% main figure, plots and title %%%%%%%%%%%%%%%%%%%%%%%  

%set(fh,'Toolbar','figure');
set(fh,'Name',stack1name);              
axes('Parent',fh,'Units','characters','Position',[10,6.15,102.4,39.4],'Tag','ha');     %main figure
        
axes('Parent',fh,'Units','characters','Position',[120 45.4 66 7.7],'Tag','xh'); 
axes('Parent',fh,'Units','characters','Position',[120 33.8 66 7.7],'Tag','yh'); 
axes('Parent',fh,'Units','characters','Position',[120 22.3 66 7.7],'Tag','zh'); 

%% %%%%%%%%%%%%%%%%%%%%%%%%% settings for the axial views %%%%%%%%%%%%%%%%
uicontrol('Units','characters',...
                'Style','text','String','profile along x',...
                'Position',[144,53.3,20,1.2]);
uicontrol('Units','characters',...
                'Style','text','String','profile along y',...
                'Position',[144,41.8,20,1.2]);
uicontrol('Units','characters',...
                'Style','text','String','profile along z',...
                'Position',[144,30.3,20,1.2]);
            
uicontrol('Units','characters',...
                'Style','text','String','number of voxels to display',...
                'Position',[188,49.2,10,3.6]); 
uicontrol('Units','characters',...
                'Style','text','String','number of voxels to display',...
                'Position',[188,37.7,10,3.6]); 
uicontrol('Units','characters',...
                'Style','text','String','number of voxels to display',...
                'Position',[188,26.2,10,3.6]); 

uicontrol('Units','characters',...
                'Tag','hxplot_range',...
                'Style','edit','String',xpix,...
                'Position',[190,47.7,6,1.2]);                  
uicontrol('Units','characters',...
                'Tag','hyplot_range',...
                'Style','edit','String',ypix,...
                'Position',[190,36.2,6,1.2]);            
uicontrol('Units','characters',...
                'Tag','hzplot_range',...
                'Style','edit','String',zpix,...
                'Position',[190,24.6,6,1.2]);

%% %%%%%%%%%%%%%%%%%%%%%%%%% zoom plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
axes('Parent',fh,'Units','characters','Position',[144 1.9 46 17.7],'Tag','zoomh');             

uicontrol('Units','characters',...
                'Style','text','String','window size',...
                'Position',[192,12,8,2.3]);
            
uicontrol('Units','characters',...
                'Tag','hzoom_width',...
                'Style','edit','String',startzoomwidth,...
                'Position',[192,10.4,6,1.2]); 

%% %%%%%%%%%%%%%%%%%%%%%%%% Threshold Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hthresh = uipanel('parent',fh,'Title','Detection Threshold',...
    'Units','characters',...
    'Position',[116,4,25,12]);

uicontrol('parent',hthresh,'Units','characters',...
            'Style','text','String','Value',...
            'Position',[1,9,8,1.2]);   
        
uicontrol('parent',hthresh,'Units','characters',...
            'Tag','hthresh_level',...
            'Style','edit','String',num2str(params.thresh.level),...
            'Position',[11,9,8,1.2]);
                
thresh = getappdata(fh,'thresh');   
if strcmp(thresh.units,'absolute')
        threshInt = thresh.level;
elseif strcmp(thresh.units,'SD')
        threshInt = thresh.level*thresh.sd;
elseif strcmp(thresh.units,'automatic')
        threshInt = thresh.thr_min + thresh.level*(thresh.thr_opt - thresh.thr_min);
end    
                    
tmp = 0;
for i=1:getappdata(fh,'nz')
    eval(['st = getappdata(fh,''st',num2str(i),''');']);
    tmp = tmp + sum(sum(st > threshInt));
end
    
setappdata(fh,'npix_above_thresh',tmp );   

uicontrol('parent',hthresh,'Units','characters',...
            'Style','text','String','# of Voxels',...
            'Position',[1,7,15,1.2]);     
        
uicontrol('parent',hthresh,'Units','characters',...
            'Tag','hnpix_above_threshold',...
            'Style','text','String',num2str(tmp),...
            'Position',[16,7,8,1.2]); 
        
if strcmp(thresh.units,'automatic')
    hth = uibuttongroup('parent',hthresh,'Units','characters',...
        'Tag','hthresh_units',...
        'Title','Units',...
        'Position',[0.5 0.2 20 5.5]);
    
        uicontrol('parent',hth,'Style','Radio','Units','characters',...
            'String','Absolute',...
            'Tag','absolute',...
            'Position',[1 3.5 18 1.5]);
        uicontrol('parent',hth,'Style','Radio','Units','characters',...
            'String','Intensity S.D.',...
            'Tag','SD',...
            'Position',[1 2 18 1.5]);
        hauto = uicontrol('parent',hth,'Style','Radio','Units','characters',...
            'String','Automatic',...
            'Tag','automatic',...
            'Position',[1 0.5 18 1.5]); 
        set(hth,'SelectedObject',hauto);
else
    hth = uibuttongroup('parent',hthresh,'Units','characters',...
        'Tag','hthresh_units',...
        'Title','Units',...
        'Position',[0.5 0.2 20 5.5]);
        hSD = uicontrol('parent',hth,'Style','Radio','Units','characters',...
            'String','Intensity S.D.',...
            'Tag','SD',...
            'Position',[1 2.5 18 1.5]);
        habs = uicontrol('parent',hth,'Style','Radio','Units','characters',...
            'String','Absolute',...
            'Tag','absolute',...
            'Position',[1 0.5 18 1.5]); 
        if strcmp(thresh.units,'SD')
             set(hth,'SelectedObject',hSD);
        else
            set(hth,'SelectedObject',habs);
        end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%% contrast panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ph = uipanel(fh,'Title','Contrast','Units','characters',...
             'Position',[30 46.2 83 8.5],'TitlePosition','centertop');

%info relative to the current slice        
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','current z-plane',...
                'Position',[10,6.5,18,1.2]);        

uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','min',...
                'Position',[0,5.4,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','median',...
                'Position',[13,5.4,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','max',...
                'Position',[24,5.4,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','SD',...
                'Position',[39,5.4,8,1.2]);   
            
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hslice_Imin_Info_val',...
                'Style','text','String',strcurImin,...
                'Position',[0,4.6,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hslice_Imed_Info_val',...
                'Style','text','String',strcurImed,...
                'Position',[13,4.6,8,1.2]);            
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hslice_Imax_Info_val',...
                'Style','text','String',strcurImax,...
                'Position',[24,4.6,8,1.2]);  
uicontrol('Parent',ph,'Units','characters',...
                'Tag','hslice_Isd_Info_val',...
                'Style','text','String',strcurIsd,...
                'Position',[39,4.6,8,1.2]);  
            
%info on whole stack                        
uicontrol('Parent',ph,'Units','characters',...
                'Style','text','String','whole stack',...
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
                'Style','pushbutton',...
                'Tag','hslice_auto',...
                'String','auto adjust on current z-plane',...
                'Units','characters','Position',[52,3,30,1.9]); 

uicontrol('Parent',ph,'Units','characters',...
                'Tag','hglob_auto',...
                'Style','pushbutton',...
                'String','auto adjust on whole stack',...
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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% z browser %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
uicontrol(fh,'Style','edit',...
                'Tag','hzpos',...
                'Units','characters',...
                'String',strZpos,...
                'Position',[54 3.5 6 1.2]);

uicontrol(fh,'Style','slider',...
                'Tag','zslider',...
                'Units','characters',...
                'Max',nz,'Min',1,'Value',1,...
                'SliderStep',[1/(nz-1) 1/(nz-1)],...
                'Position',[10 2.3 102.4 1.2]);
%% action button
uicontrol(fh,'Style','pushbutton',...
                    'Tag','hclose',...
                    'Units','characters',...
                    'String','done',...
                    'Position',[116 0.4 22 3]);
%% %%%%%%%%%%%%%%%%%%% overlay %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

uicontrol(fh,'Style','togglebutton',...
                    'Tag','hoverlay',...
                    'Units','characters',...
                    'Value',0,...
                    'String','Hide Overlay',...
                    'Position',[50 0.2 22 1.5]);                
                
                             
clear('nz','strXpos','strYpos','strZpos',...
    'Imed','Isd','Imin','Imax',...
    'strImin','strImax','strImed','strIsd',...
    'curImin','curImax','curImed','curIsd',...
    'strcurImin','strcurImax','strcurImed','strcurIsd',...
    'xpix','ypix','zpix',...
    'startzoomwidth','ph');

clear('stack1','stack1name');
end

