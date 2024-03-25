function pos2 = view_stack_col2(varargin)
% VIEW_STACK If you need help call 911 .
% Timothee lionnet december 6, 2008

%one argument: one color channel (data)
%two arguments: the first one is the data (green channel)
%the second one is the overlay (red channel) - typically the position of
%the spots located by some detection algorithm
%three arguments: the first one is the data (green channel)
%the second one is the overlay (red channel) - typically the position of
%the spots located by some detection algorithm
%the third one is a string specifying the shape of the object locating the spots in the overlay channel

%position of x and y of the pointer when I freeze the xy browsing
global nx ny nz stack2 spots_range pos2

%% optional arguments that specify 
%1) the shape of the red overlay marker
dotshape = 'square';
%2)how to treat the overlay (sparse spots vs. normal overlay)
overlay_mode = 'simple';    %simple overlay.

%overlay_mode = 'spots';    %allows user to select spots depending on their
%values and positions. Not recommended for arrays with large numbers of non
%zeros values in the red channel

if nargin >=1
    for i = 1:nargin
        if ischar(varargin{i})
            if strcmp(varargin{i},'dotshape')
                if nargin >= i+1
                    dotshape = varargin{i+1};
                end
            elseif strcmp(varargin{i},'overlay_mode')
                if nargin >= i+1
                    overlay_mode = varargin{i+1};
                    dotshape = 'small';
                end
            end
        end

    end
end

%% the main data (green channel)
if nargin == 0
    stack1 = timtiffread();
elseif nargin >= 1
    stack1 = varargin{1};
end

stack1 = double(stack1); 
if ndims(stack1) ~= 3 || size(stack1,3)==1
    fprintf('I like my arguments to be 3 dimensional, babe.\n');
    return;
end
[nx ny nz] = size(stack1);
stack1name = inputname(1);

%% checking whether we have two 3D stacks - this is the overlay (red channel)
if nargin>=2
    if isnumeric(varargin{2})
        if ndims(varargin{2}) == 2 && size(varargin{2},2) ==4
            pos = varargin{2};
            stack2  = convert_coordinates_to_stack(pos,nx,ny,nz);
        elseif ndims(varargin{2}) == 3 
            stack2 = varargin{2};
            if size(stack1) ~= size(stack2)
                fprintf('I like my stacks to be the same size, babe.\n');
                return;
            end
            if strcmp(overlay_mode,'spots')
                pos = convert_stack_to_coordinates(stack2);
            else
                pos = 0;
            end
        else
            fprintf('I like my arguments to have the right size, babe.\n');
            return;
        end
    end
else
    stack2 = zeros(nx,ny,nz);
    pos = 0;
end
   
pos2 = pos;

if strcmp('small',dotshape)
    %dots appear as they are in stack2
elseif strcmp('large',dotshape)
    stack2 = dilate(stack2,'cube'); %spots are expanded as 3x3x3 cubes to help visualization
elseif strcmp('mix',dotshape)
    stack2 = dilate(stack2,'mix');  %spots are expanded as 3x3 square + cross to help visualization
elseif strcmp('square',dotshape)
    stack2 = dilate(stack2,'square'); %spots are expanded as 3x3 square
end 

stack2val = stack2; %keeps the values of the second array so that they are tracable on the display panel
stack2 = double(stack2);
max2 = max(max(max(stack2))); min2 = min(min(min(stack2)));
if max2 == min2
    stack2 = 0*stack2;
else
    stack2(:,:,:) = max(  min( (stack2(:,:,:) - min2 )/(max2-min2) , 1 ) , 0) ; %overlay data should be between 0 and 1. Convolution can lead to >1 values.
end
%% setting the different panels
[fh,ha,penseedujour,hoverlay,hclose,...
    xh,yh,zh,hxplot_title,hyplot_title,hzplot_title,hxplot_size_txt,hyplot_size_txt,hzplot_size_txt,...
    hxplot_range,hyplot_range,hzplot_range,zoomh,hzoom_title,hzoom_width,...
    ph,hslice_title,hslice_Imin_Info_title,hslice_Imed_Info_title,hslice_Imax_Info_title,hslice_Isd_Info_title,...
    hslice_Imin_Info_val,hslice_Imed_Info_val,hslice_Imax_Info_val,hslice_Isd_Info_val,...
    hglobal_Imin_Info_title,hglobal_Imed_Info_title,hglobal_Imax_Info_title,hglobal_Isd_Info_title,...
    hglobal_Imin_Info_val,hglobal_Imed_Info_val,hglobal_Imax_Info_val,hglobal_Isd_Info_val,...
    hslice_auto,hglob_auto,hImintitle,hImin,hImaxtitle,hImax,...
    hstate,hxtitle,hytitle,hxpos,hypos,hI1title,hI2title,hI1val,hI2val,hzpos,zslider,...
    spots_range,...
    hspots_title,hnumber_of_spots_selected,hnumber_of_spots_selected_title,...
    htotal_number_of_spots,htotal_number_of_spots_title,...
    hspots_xrange_title,hspots_yrange_title,hspots_zrange_title,hspots_Irange_title,...
    hspots_xmin,hspots_ymin,hspots_zmin,hspots_Imin,...
    hspots_xmax,hspots_ymax,hspots_zmax,hspots_Imax,...
    hgauss_title,h_x0_title,h_x0,h_y0_title,h_y0,h_z0_title,h_z0,...
    h_I0_title,h_I0,h_bg0_title,h_bg0,...
    h_sxy0_title,h_sxy0,h_sz0_title,h_sz0] = set_figure_and_panel_for_view_stack_col2(stack1,stack2,pos,overlay_mode,stack1name);

%% setting all the callback functions       
set(fh,'WindowButtonMotionFcn', {@viewlocaldata,stack1,stack2val,hxpos,hypos,...
                    hI1val,hI2val,zslider,hImin,hImax,xh,hxplot_range,...
                    yh,hyplot_range,zh,hzplot_range,hzoom_width,zoomh,ha,hstate,hoverlay});
set(fh,'WindowButtonUpFcn', {@change_state,hstate,ha});   
set(hclose,'Callback', {@clean_close,fh});   
set(zslider,'Callback',{@slider_change,hzpos,hI1val,hI2val,stack1,stack2val,hImin,hImax,...
            hslice_Imin_Info_val,hslice_Imed_Info_val,hslice_Imax_Info_val,hslice_Isd_Info_val,ha,fh,hoverlay,...
            xh,yh,zh,zoomh,hxplot_range,hyplot_range,hzplot_range,hzoom_width,hstate});
set(hzpos,'Callback',{@zpos_change,nz,zslider,hI1val,hI2val,stack1,stack2val,hImin,hImax,...
            hslice_Imin_Info_val,hslice_Imed_Info_val,hslice_Imax_Info_val,hslice_Isd_Info_val,ha,fh,hoverlay,...
            xh,yh,zh,zoomh,hxplot_range,hyplot_range,hzplot_range,hzoom_width,hstate});            
set(hImin,'Callback',{@change_Iminmax,zslider,stack1,hImin,hImax,ha,fh,hoverlay,hzoom_width,hstate,zoomh});
set(hImax,'Callback',{@change_Iminmax,zslider,stack1,hImin,hImax,ha,fh,hoverlay,hzoom_width,hstate,zoomh});
set(hslice_auto,'Callback',{@slice_auto_adjust,zslider,stack1,hImin,hImax,ha,fh,hoverlay,hzoom_width,hstate,zoomh});
set(hglob_auto,'Callback',{@global_auto_adjust,zslider,stack1,hImin,hImax,ha,fh,hoverlay,hzoom_width,hstate,zoomh});
set(hoverlay,'Callback',{@switch_overlay,ha,fh,zoomh,hImin,hImax,hzpos,stack1,hzoom_width,hstate});

if strcmp(overlay_mode,'spots')
    set(hspots_xmin,'Callback',{@change_spots_range,ha,fh,zoomh,hImin,hImax,zslider,stack1,pos,hzoom_width,hstate,hoverlay,hnumber_of_spots_selected,'xmin'});
    set(hspots_xmax,'Callback',{@change_spots_range,ha,fh,zoomh,hImin,hImax,zslider,stack1,pos,hzoom_width,hstate,hoverlay,hnumber_of_spots_selected,'xmax'});
    set(hspots_ymin,'Callback',{@change_spots_range,ha,fh,zoomh,hImin,hImax,zslider,stack1,pos,hzoom_width,hstate,hoverlay,hnumber_of_spots_selected,'ymin'});
    set(hspots_ymax,'Callback',{@change_spots_range,ha,fh,zoomh,hImin,hImax,zslider,stack1,pos,hzoom_width,hstate,hoverlay,hnumber_of_spots_selected,'ymax'});
    set(hspots_zmin,'Callback',{@change_spots_range,ha,fh,zoomh,hImin,hImax,zslider,stack1,pos,hzoom_width,hstate,hoverlay,hnumber_of_spots_selected,'zmin'});
    set(hspots_zmax,'Callback',{@change_spots_range,ha,fh,zoomh,hImin,hImax,zslider,stack1,pos,hzoom_width,hstate,hoverlay,hnumber_of_spots_selected,'zmax'});
    set(hspots_Imin,'Callback',{@change_spots_range,ha,fh,zoomh,hImin,hImax,zslider,stack1,pos,hzoom_width,hstate,hoverlay,hnumber_of_spots_selected,'Imin'});
    set(hspots_Imax,'Callback',{@change_spots_range,ha,fh,zoomh,hImin,hImax,zslider,stack1,pos,hzoom_width,hstate,hoverlay,hnumber_of_spots_selected,'Imax'});

elseif strcmp(overlay_mode,'simple')
    set(hgauss_title,'Callback',{@local_3D_gaussian_fit_from_viewer,hzpos,stack1,h_x0,h_y0,h_z0,h_I0,h_bg0,h_sxy0,h_sz0});
end

set(fh,'SelectionType','alt');

plot_current_z_stack(fh,ha,hImin,hImax,stack1,stack2,1,hoverlay);
set(fh,'Visible','on');
clear('temp','nx','ny','nz');
uiwait;
clear('stack1','stack2','stack2val','max2','min2','pos',...
    'dotshape','overlay_mode','stack1name');
clear('fh','ha','penseedujour','hoverlay','hclose',...
    'xh','yh','zh','hxplot_title','hyplot_title','hzplot_title','hxplot_size_txt','hyplot_size_txt','hzplot_size_txt',...
    'hxplot_range','hyplot_range','hzplot_range','zoomh','hzoom_title','hzoom_width',...
    'ph','hslice_title','hslice_Imin_Info_title','hslice_Imed_Info_title','hslice_Imax_Info_title','hslice_Isd_Info_title',...
    'hslice_Imin_Info_val','hslice_Imed_Info_val','hslice_Imax_Info_val','hslice_Isd_Info_val',...
    'hglobal_Imin_Info_title','hglobal_Imed_Info_title','hglobal_Imax_Info_title','hglobal_Isd_Info_title',...
    'hglobal_Imin_Info_val','hglobal_Imed_Info_val','hglobal_Imax_Info_val','hglobal_Isd_Info_val',...
    'hslice_auto','hglob_auto','hImintitle','hImin','hImaxtitle','hImax',...
    'hstate','hxtitle','hytitle','hxpos','hypos','hI1title','hI2title','hI1val','hI2val','hzpos','zslider',...
    'spots_range',...
    'hspots_title','hnumber_of_spots_selected','hnumber_of_spots_selected_title',...
    'htotal_number_of_spots','htotal_number_of_spots_title',...
    'hspots_xrange_title','hspots_yrange_title','hspots_zrange_title','hspots_Irange_title',...
    'hspots_xmin','hspots_ymin','hspots_zmin','hspots_Imin',...
    'hspots_xmax','hspots_ymax','hspots_zmax','hspots_Imax',...
    'hgauss_title','h_x0_title','h_x0','h_y0_title','h_y0','h_z0_title','h_z0',...
    'h_I0_title','h_I0','h_bg0_title','h_bg0',...
    'h_sxy0_title','h_sxy0','h_sz0_title','h_sz0');
return;

end

%%%%%%%%%%%%%%%%%%%%%%%%%% Callback functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function local_3D_gaussian_fit_from_viewer(src,eventdata,hzpos,stack1,h_x0,h_y0,h_z0,h_I0,h_bg0,h_sxy0,h_sz0)

    global xfreeze yfreeze
    p = 0;
    
    z = str2num(get(hzpos,'String'));
    z = round(z);
    
    [Gaussout, resnorm]=gaussian_3D_fit_local(stack1,xfreeze,yfreeze,z,p);
    
    set(h_I0,'String',num2str(Gaussout(1)));
    set(h_bg0,'String',num2str(Gaussout(2)));
    set(h_sxy0,'String',num2str(Gaussout(3)));
    set(h_sz0,'String',num2str(Gaussout(4)));
    set(h_x0,'String',num2str(Gaussout(5)));
    set(h_y0,'String',num2str(Gaussout(6)));
    set(h_z0,'String',num2str(Gaussout(7)));
    
    clear('xfreeze','yfreeze','z','p');
    clear('src','eventdata','hzpos','stack1','h_x0','h_y0','h_z0','h_I0','h_bg0','h_sxy0','h_sz0');
end

function change_spots_range(src,eventdata,ha,fh,zoomh,hImin,hImax,zslider,stack1,pos,hzoom_width,hstate,hoverlay,hnumber_of_spots_selected,mode)
    
    global stack2 spots_range pos2
    [nx ny nz] = size(stack1);
    
    if strcmp(mode(1),'x'), ncol =1; end
    if strcmp(mode(1),'y'), ncol =2; end
    if strcmp(mode(1),'z'), ncol =3; end
    if strcmp(mode(1),'I'), ncol =4; end
         
    if strcmp(mode(2:4),'min')
            spots_range(1,ncol) = str2double(get(src,'String'));
    elseif strcmp(mode(2:4),'max')  
            spots_range(2,ncol) = str2double(get(src,'String'));
    end
    pos2 = select_array_by_multiple_col_values(pos,spots_range);
    if pos2 == 0
        nspots = 0;
    else
        nspots = num2str(size(pos2,1));
    end
    set(hnumber_of_spots_selected,'String',nspots);
    stack2  = convert_coordinates_to_stack(pos2,nx,ny,nz,'ones');
    
    z = get(zslider,'Value');
    %replotting the main window
    plot_current_z_stack(fh,ha,hImin,hImax,stack1,stack2,z,hoverlay);
    
    %replotting the zoom window
    zoomwidth = str2double(get(hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh,ha,stack1,hstate); %I get the current (x,y) of the mouse
    
    plot_zoom(x,y,z,zoomwidth,zoomh,stack1,stack2,fh,hImin,hImax,hoverlay);
    clear('x','y','z','nx','ny','nz','zoomwidth','Imax','Imin','ncol','nspots');
    clear('stack2','spots_range','pos2');
    clear('src','eventdata','ha','fh','zoomh','hImin','hImax','zslider','stack1',...
        'pos','hzoom_width','hstate','hoverlay','hnumber_of_spots_selected','mode');
end

function clean_close(src,eventdata,fh)
close(fh);
clear('src','eventdata','fh');
end

function change_state(src,eventdata,hstate,ha)
    %switches between the snap and grab mode
    %hstate value = 1: grap
    %hstate value = 2: snap
    global xfreeze yfreeze
    
    val = get(hstate,'Value');
    val = 3-val; %1->2 or 2->1
    set(hstate,'Value',val);
    
    if val == 2   %if switching to snap state
        set(src,'CurrentAxes',ha);
        point = get(gca,'CurrentPoint');
        yfreeze = round(point(1,1));    %note the inverse convention for mouse and array (xy -> yx)
        xfreeze = round(point(1,2));
    end
    clear('point','val','xfreeze','yfreeze');
    clear('src','eventdata','hstate','ha');
end

function switch_overlay(hoverlay,eventdata,ha,fh,zoomh,hImin,hImax,hzpos,stack1,hzoom_width,hstate)
    global xfreeze yfreeze stack2
    [nx ny nz] = size(stack1);
    z = str2num(get(hzpos,'String'));
    z = round(z);
    
    %switching the state of the button
    is_overlay_on = 1 - get(hoverlay,'Value');
    if is_overlay_on == 1 
        set(hoverlay,'String','hide overlay');
    else 
        set(hoverlay,'String','show overlay');
    end
    
    %replotting the zoom window
    zoomwidth = str2num(get(hzoom_width,'String'));
    
    val = get(hstate,'Value');  %this is the snap/grab mode for the local data panel
    if val == 1                 %if grab mode is activated, I replot the zoom data at some arbitrary location
        plot_zoom(1,ny,z,zoomwidth,zoomh,stack1,stack2,fh,hImin,hImax,hoverlay);
          
    else                %if snap mode is activated, I plot the local data at the same location
        plot_zoom(xfreeze,yfreeze,z,zoomwidth,zoomh,stack1,stack2,fh,hImin,hImax,hoverlay);
    end
    
    %replotting the main window
    plot_current_z_stack(fh,ha,hImin,hImax,stack1,stack2,z,hoverlay);
    
    clear('z','val','zoomwidth','is_overlay_on');
    clear('nx','ny','nz','xfreeze','yfreeze','stack2');
    clear('hoverlay','eventdata','ha','fh','zoomh','hImin','hImax','hzpos','stack1','hzoom_width','hstate');
end

function slider_change(zslider,eventdata,hzpos,hI1val,hI2val,stack1,stack2val,hImin,hImax,...
    hslice_Imin_Info_val,hslice_Imed_Info_val,hslice_Imax_Info_val,hslice_Isd_Info_val,ha,fh,hoverlay,...
    xh,yh,zh,zoomh,hxplot_range,hyplot_range,hzplot_range,hzoom_width,hstate)
    
    global xfreeze yfreeze stack2
    
    [nx ny nz] = size(stack1);
    %updating the z-slider and the z info window
    z = get(zslider,'Value');
    z = round(z);
    set(zslider,'Value',z);
    set(hzpos,'String',num2str(z));
    
    %updating the Imin/Imax info
    temp = reshape(stack1(1:nx,1:ny,z),1,nx*ny);
    curImin = min(temp); curImax = max(temp); curImed = median(temp); curIsd = std(temp);
    clear('temp');
    set(hslice_Imin_Info_val,'String',num2str(curImin));
    set(hslice_Imed_Info_val,'String',num2str(curImed));
    set(hslice_Imax_Info_val,'String',num2str(curImax)); 
    set(hslice_Isd_Info_val,'String',num2str(curIsd)); 
    
    %replotting the local data (zoom + axial views)
    [x,y] = get_current_pointer_position(fh,ha,stack1,hstate); %I get the current (x,y) of the mouse
    Int1 = stack1(x,y,z);
    Int2 = stack2val(x,y,z);
    
    set(hI1val,'String',num2str(Int1));
    set(hI2val,'String',num2str(Int2));

    plot_local_data(x,y,z,fh,hImin,hImax,xh,yh,zh,zoomh,...
        stack1,hxplot_range,hyplot_range,hzplot_range,hzoom_width,hoverlay);
    
    %replotting the main window
    plot_current_z_stack(fh,ha,hImin,hImax,stack1,stack2,z,hoverlay);  
    
    clear('x','y','z','temp','curImin','curImax','curImed');
    clear('nx','ny','nz','xfreeze','yfreeze','stack2');
    clear('Int1','Int2');
    clear('zslider','eventdata','hzpos','hI1val','hI2val','stack1','stack2val','hImin','hImax',...
    'hslice_Imin_Info_val','hslice_Imed_Info_val','hslice_Imax_Info_val','hslice_Isd_Info_val','ha','fh','hoverlay',...
    'xh','yh','zh','zoomh','hxplot_range','hyplot_range','hzplot_range','hzoom_width','hstate');
end

function zpos_change(hzpos,eventdata,nz,zslider,hI1val,hI2val,stack1,stack2val,hImin,hImax,...
    hslice_Imin_Info_val,hslice_Imed_Info_val,hslice_Imax_Info_val,hslice_Isd_Info_val,ha,fh,hoverlay,...
    xh,yh,zh,zoomh,hxplot_range,hyplot_range,hzplot_range,hzoom_width,hstate)

    global xfreeze yfreeze stack2
    
    [nx ny nz] = size(stack1);
    %updating the z-slider and the z info window
    z = str2num(get(hzpos,'String'));
    z = round(z);
    z = max(min(z,nz) , 1);
    set(hzpos,'String',num2str(z));
    set(zslider,'Value',z);
    
    %updating the Imin/Imax info
    temp = reshape(stack1(1:nx,1:ny,2),1,nx*ny);
    curImin = min(temp); curImax = max(temp); curImed = median(temp); curIsd = std(temp);
    clear('temp');
    set(hslice_Imin_Info_val,'String',num2str(curImin));
    set(hslice_Imed_Info_val,'String',num2str(curImed));
    set(hslice_Imax_Info_val,'String',num2str(curImax));
    set(hslice_Isd_Info_val,'String',num2str(curIsd)); 
    
    %replotting the main window
    plot_current_z_stack(fh,ha,hImin,hImax,stack1,stack2,z,hoverlay);
    
    %replotting the local data (zoom + axial views)
    [x,y] = get_current_pointer_position(fh,ha,stack1,hstate); %I get the current (x,y) of the mouse
    Int1 = stack1(x,y,z);
    set(hI1val,'String',num2str(Int1));
    Int2 = stack2val(x,y,z);
    set(hI2val,'String',num2str(Int2));

    plot_local_data(x,y,z,fh,hImin,hImax,xh,yh,zh,zoomh,...
        stack1,hxplot_range,hyplot_range,hzplot_range,hzoom_width,hoverlay);
    
    clear('x','y','z','curImin','curImax','curImed','curIsd'); 
    clear('nx','ny','nz','xfreeze','yfreeze','stack2','Int1','Int2');
    clear('hzpos','eventdata','nz','zslider','hI1val','hI2val','stack1','stack2val','hImin','hImax',...
    'hslice_Imin_Info_val','hslice_Imed_Info_val','hslice_Imax_Info_val','hslice_Isd_Info_val','ha','fh','hoverlay',...
    'xh','yh','zh','zoomh','hxplot_range','hyplot_range','hzplot_range','hzoom_width','hstate');
end

function change_Iminmax(hObj,eventdata,zslider,stack1,hImin,hImax,ha,fh,hoverlay,hzoom_width,hstate,zoomh)
    
    global stack2
    %updating the Imin/Imax info
    Imin = str2num(get(hImin,'String'));
    Imax = str2num(get(hImax,'String'));
    
    z = get(zslider,'Value');
    
    %replotting the main window
    plot_current_z_stack(fh,ha,hImin,hImax,stack1,stack2,z,hoverlay);
    
    %replotting the zoom window
    zoomwidth = str2num(get(hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh,ha,stack1,hstate); %I get the current (x,y) of the mouse
    
    plot_zoom(x,y,z,zoomwidth,zoomh,stack1,stack2,fh,hImin,hImax,hoverlay);
    clear('x','y','z','zoomwidth','Imax','Imin','stack2','zoomwidth');
    clear('hObj','eventdata','zslider','stack1','hImin','hImax','ha','fh','hoverlay','hzoom_width','hstate','zoomh');
end

function slice_auto_adjust(hslice_auto,eventdata,zslider,stack1,hImin,hImax,ha,fh,hoverlay,hzoom_width,hstate,zoomh)
    global stack2
    z = get(zslider,'Value');
    
    %updating Imin / Imax info
    Imin = min(min(stack1(:,:,z))); Imax = max(max(stack1(:,:,z)));
    set(hImin,'String',num2str(Imin)); set(hImax,'String',num2str(Imax));
    
    %replotting the main window
    plot_current_z_stack(fh,ha,hImin,hImax,stack1,stack2,z,hoverlay);
    
    %replotting the zoom window
    zoomwidth = str2num(get(hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh,ha,stack1,hstate); %I get the current (x,y) of the mouse
    
    plot_zoom(x,y,z,zoomwidth,zoomh,stack1,stack2,fh,hImin,hImax,hoverlay);
    clear('x','y','z','zoomwidth','Imin','Imax','stack2','zoomwidth');
    clear('hslice_auto','eventdata','zslider','stack1','hImin','hImax','ha','fh','hoverlay','hzoom_width','hstate','zoomh');
end

function global_auto_adjust(hglob_auto,eventdata,zslider,stack1,hImin,hImax,ha,fh,hoverlay,hzoom_width,hstate,zoomh)
    global stack2
    z = get(zslider,'Value');
    
    %updating the Imin/Imax info
    Imin = min(min(min(stack1))); Imax = max(max(max(stack1)));
    set(hImin,'String',num2str(Imin)); set(hImax,'String',num2str(Imax));
    
    %replotting the main window
    plot_current_z_stack(fh,ha,hImin,hImax,stack1,stack2,z,hoverlay);
    
    %replotting the zoom window
    zoomwidth = str2num(get(hzoom_width,'String'));
    [x,y] = get_current_pointer_position(fh,ha,stack1,hstate); %I get the current (x,y) of the mouse
    
    plot_zoom(x,y,z,zoomwidth,zoomh,stack1,stack2,fh,hImin,hImax,hoverlay);
    
    clear('x','y','z','zoomwidth','Imin','Imax','stack2');
    clear('hglob_auto','eventdata','zslider','stack1','hImin','hImax','ha','fh','hoverlay','hzoom_width','hstate','zoomh');
end

function viewlocaldata(src,eventdata,stack1,stack2val,hxpos,hypos,hI1val,hI2val,zslider,hImin,hImax,...
    xh,hxplot_range,yh,hyplot_range,zh,hzplot_range,hzoom_width,zoomh,ha,hstate,hoverlay)
    
    global stack2
    [nx ny nz] = size(stack1);
    val = get(hstate,'Value');          %this is the snap/grab mode
    
    %I replot things only if the grab mode is active
    if val == 1 

        %I get the current (x,y,z) 
        [x,y] = get_current_pointer_position(src,ha,stack1,hstate); 
        z = get(zslider,'Value');
        z = round(z);

        %I update the current display of the pointer info
        set(hxpos,'String',num2str(x));
        set(hypos,'String',num2str(y));

        if(x>=1 && y >=1 && x <= nx && y <= ny)
            Int1 = stack1(x,y,z);
            set(hI1val,'String',num2str(Int1));
            Int2 = stack2val(x,y,z);
            set(hI2val,'String',num2str(Int2));

            %plotting the profiles and zoom window
            plot_local_data(x,y,z,src,hImin,hImax,xh,yh,zh,zoomh,...
            stack1,hxplot_range,hyplot_range,hzplot_range,hzoom_width,hoverlay);
        end
    end
    clear('x','y','z','val','Int1','Int2','nx','ny','nz','stack2');
    clear('src','eventdata','stack1','stack2val','hxpos','hypos','hI1val','hI2val','zslider','hImin','hImax',...
    'xh','hxplot_range','yh','hyplot_range','zh','hzplot_range','hzoom_width','zoomh','ha','hstate','hoverlay');
end

function plot_local_data(x,y,z,fh,hImin,hImax,xh,yh,zh,zoomh,...
    stack1,hxplot_range,hyplot_range,hzplot_range,hzoom_width,hoverlay)
    
    global stack2
    xrange = str2num(get(hxplot_range,'String'));
    yrange = str2num(get(hyplot_range,'String'));
    zrange = str2num(get(hzplot_range,'String'));
    zoomwidth = str2num(get(hzoom_width,'String'));

    %plotting the profiles and zoom window
    plot_xprofile(x,y,z,xrange,stack1,xh,fh);
    plot_yprofile(x,y,z,yrange,stack1,yh,fh);
    plot_zprofile(x,y,z,zrange,stack1,zh,fh);
    plot_zoom(x,y,z,zoomwidth,zoomh,stack1,stack2,fh,hImin,hImax,hoverlay);
    
    clear('xrange','yrange','zrange','zoomwidth','stack2');
    clear('x','y','z','fh','hImin','hImax','xh','yh','zh','zoomh',...
    'stack1','hxplot_range','hyplot_range','hzplot_range','hzoom_width','hoverlay');
end

function [xarr yarr] = get_current_pointer_position(fh,ha,stack,hstate)

global xfreeze yfreeze

[nx ny nz] = size(stack);
val = get(hstate,'Value');%this is the snap/grab mode for the local data panel

if val == 1         %if grab mode is activated, I replot the local data at some arbitrary location
        set(fh,'CurrentAxes',ha);
        point = get(gca,'CurrentPoint');
        xpic = round(point(1,1));   ypic = round(point(1,2));
        xarr = ypic;    yarr = xpic;    %pic and array have different xy conventions
        xarr = max(xarr,1);  xarr = min(xarr,nx);
        yarr = max(yarr,1);  yarr = min(yarr,ny);
else
        xarr = xfreeze; yarr = yfreeze;
end

clear('nx','ny','nz','val','point','xpic','ypic','xfreeze','yfreeze');
clear('fh','ha','stack','hstate');
end

function plot_xprofile(xarr,yarr,z,xrange,stack1,xh,fh)
nx = size(stack1,1);
xprofile(1:nx) = stack1(1:nx,yarr,z);

if 2*xrange+1>nx
        set(fh,'CurrentAxes',xh);
        Imin = min(xprofile(1:nx));
        Imax = max(xprofile(1:nx));
        plot(1:nx,xprofile(1:nx),[xarr xarr],[Imin Imax]); 
elseif (xarr - xrange >= 1) && (xarr + xrange <= nx)
        set(fh,'CurrentAxes',xh);
        Imin = min(xprofile(xarr-xrange:xarr+xrange));
        Imax = max(xprofile(xarr-xrange:xarr+xrange));
        plot(xarr-xrange:xarr+xrange,xprofile(xarr-xrange:xarr+xrange),...
            [xarr xarr],[Imin Imax]);
elseif xarr - xrange < 1 && xarr + xrange <= nx 
        set(fh,'CurrentAxes',xh);
        x = xrange+1;
        Imin = min(xprofile(x-xrange:x+xrange));
        Imax = max(xprofile(x-xrange:x+xrange));
        plot(x-xrange:x+xrange,xprofile(x-xrange:x+xrange),...
            [xarr xarr],[Imin Imax]);
elseif xarr + xrange > nx && xarr - xrange >= 1
        set(fh,'CurrentAxes',xh);
        x = nx-xrange;
        Imin = min(xprofile(x-xrange:x+xrange));
        Imax = max(xprofile(x-xrange:x+xrange));
        plot(x-xrange:x+xrange,xprofile(x-xrange:x+xrange),...
            [xarr xarr],[Imin Imax]);   
end 
clear('nx','xprofile','Imin','Imax','x');
clear('xarr','yarr','z','xrange','stack1','xh','fh');
end

function plot_yprofile(xarr,yarr,z,yrange,stack1,yh,fh)
ny = size(stack1,2);
yprofile(1:ny) = stack1(xarr,1:ny,z);

if 2*yrange+1>ny
        set(fh,'CurrentAxes',yh);
        Imin = min(yprofile(1:ny));
        Imax = max(yprofile(1:ny));
        plot(1:ny,yprofile(1:ny),[yarr yarr],[Imin Imax]);
elseif yarr - yrange >= 1 && yarr + yrange <= ny
        set(fh,'CurrentAxes',yh);
        Imin = min(yprofile(yarr-yrange:yarr+yrange));
        Imax = max(yprofile(yarr-yrange:yarr+yrange));
        plot(yarr-yrange:yarr+yrange,yprofile(yarr-yrange:yarr+yrange),...
            [yarr yarr],[Imin Imax]);
elseif yarr - yrange < 1  && yarr + yrange <= ny
        set(fh,'CurrentAxes',yh);
        y = yrange+1;
        Imin = min(yprofile(y-yrange:y+yrange));
        Imax = max(yprofile(y-yrange:y+yrange));
        plot(y-yrange:y+yrange,yprofile(y-yrange:y+yrange),...
            [yarr yarr],[Imin Imax]);
elseif yarr + yrange > ny && yarr - yrange >= 1
        set(fh,'CurrentAxes',yh);
        y = ny-yrange;
        Imin = min(yprofile(y-yrange:y+yrange));
        Imax = max(yprofile(y-yrange:y+yrange));
        plot(y-yrange:y+yrange,yprofile(y-yrange:y+yrange),...
            [yarr yarr],[Imin Imax]);
end 
clear('ny','yprofile','Imin','Imax','y');
clear('xarr','yarr','z','yrange','stack1','yh','fh');
end

function plot_zprofile(xarr,yarr,z,zrange,stack1,zh,fh)
nz = size(stack1,3);
zprofile(1:nz) = stack1(xarr,yarr,1:nz);

if 2*zrange+1 > nz
        set(fh,'CurrentAxes',zh);
        Imin = min(zprofile(1:nz));
        Imax = max(zprofile(1:nz));
        plot(1:nz,zprofile(1:nz),[z z],[Imin Imax]);
elseif z - zrange >= 1 && z + zrange <= nz
        set(fh,'CurrentAxes',zh);
        Imin = min(zprofile(z-zrange:z+zrange));
        Imax = max(zprofile(z-zrange:z+zrange));
        plot(z-zrange:z+zrange,zprofile(z-zrange:z+zrange),[z z],[Imin Imax]);
elseif z - zrange < 1 && z + zrange <= nz  
        set(fh,'CurrentAxes',zh);
        zz = zrange+1;
        Imin = min(zprofile(zz-zrange:zz+zrange));
        Imax = max(zprofile(zz-zrange:zz+zrange));
        plot(zz-zrange:zz+zrange,zprofile(zz-zrange:zz+zrange),...
            [z z],[Imin Imax]);
elseif z + zrange > nz && z - zrange >= 1
        set(fh,'CurrentAxes',zh);
        zz = nz-zrange;
        Imin = min(zprofile(zz-zrange:zz+zrange));
        Imax = max(zprofile(zz-zrange:zz+zrange));
        plot(zz-zrange:zz+zrange,zprofile(zz-zrange:zz+zrange),...
            [z z],[Imin Imax]);
end 

clear('nz','zprofile','Imin','Imax','z');
clear('xarr','yarr','z','zrange','stack1','zh','fh');
end
