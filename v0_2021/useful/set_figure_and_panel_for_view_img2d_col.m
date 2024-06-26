function [fh,ha,penseedujour,hoverlay,hclose,...
    xh,yh,hxplot_title,hyplot_title,hxplot_size_txt,hyplot_size_txt,...
    hxplot_range,hyplot_range,zoomh,hzoom_title,hzoom_width,...
    ph,...
    hglobal_Imin_Info_title,hglobal_Imed_Info_title,hglobal_Imax_Info_title,hglobal_Isd_Info_title,...
    hglobal_Imin_Info_val,hglobal_Imed_Info_val,hglobal_Imax_Info_val,hglobal_Isd_Info_val,...
    hglob_auto,hImintitle,hImin,hImaxtitle,hImax,...
    hstate,hxtitle,hytitle,hxpos,hypos,hI1title,hI2title,hI1val,hI2val,...
    spots_range,...
    hspots_title,hnumber_of_spots_selected,hnumber_of_spots_selected_title,...
    htotal_number_of_spots,htotal_number_of_spots_title,...
    hspots_xrange_title,hspots_yrange_title,hspots_Irange_title,...
    hspots_xmin,hspots_ymin,hspots_Imin,...
    hspots_xmax,hspots_ymax,hspots_Imax,...
    hgauss_title,h_x0_title,h_x0,h_y0_title,h_y0,...
    h_I0_title,h_I0,h_bg0_title,h_bg0,...
    h_sxy0_title,h_sxy0] = set_figure_and_panel_for_view_img2d_col(stack1,stack2,pos,overlay_mode,stack1name)

%% main figure, plots and title 
fh = figure(...
              'Units','characters',...
              'MenuBar','none',...
              'Toolbar','none',...
              'Position',[20 3 200 55],...
              'Visible','off');   

set(fh,'Toolbar','figure');
set(fh,'Name',stack1name);  

ha = axes('Units','Pixels','Position',[50,80,512,512]);     %main figure

penseedujour = uicontrol('Style','text','String','another day another dollar',...
           'Position',[50,670,80,30]);
        
hoverlay = uicontrol('Units','pixels',...
                'Style','togglebutton',...
                'String','hide overlay',...
                'Units','pixels','Position',[200,5,170,15]);

hclose = uicontrol(fh,'Style','pushbutton',...
                    'String','close',...
                    'Position',[950 5 50,15]);              
  
%% initialize variables            
[nx,ny,strXpos, strYpos, strI1val, strI2val, ...
    Imin,Imax,Imed,Isd,strImin,strImax,strImed,strIsd,...
    xpix,ypix,startzoomwidth,...
    spots_Imin,spots_Imax,spots_range] = init_variables_for_view_stack_col2_img(stack1,stack2,pos,overlay_mode);          
                      
%% settings for the axial views 
[xh,yh,hxplot_title,hyplot_title,...
    hxplot_size_txt,hyplot_size_txt,...
    hxplot_range,hyplot_range] = build_axial_view_plots_and_settings_img(fh,xpix,ypix);

%% zoom plot           
[zoomh,hzoom_title,hzoom_width] = get_zoom_window_img(fh,startzoomwidth);
            
%% selection of overlay spots
% appears only when overlay_mode = 'spots'
[hspots_title,hnumber_of_spots_selected,hnumber_of_spots_selected_title,...
    htotal_number_of_spots,htotal_number_of_spots_title,...
    hspots_xrange_title,hspots_yrange_title,hspots_Irange_title,...
    hspots_xmin,hspots_ymin,hspots_Imin,...
    hspots_xmax,hspots_ymax,hspots_Imax] = get_overlay_spots_settings_img(fh,pos,spots_range,overlay_mode);

%% gaussian fit panel 
% appears only when overlay_mode = 'simple'
[hgauss_title,h_x0_title,h_x0,h_y0_title,h_y0,...
    h_I0_title,h_I0,h_bg0_title,h_bg0,...
    h_sxy0_title,h_sxy0] = get_local_gaussian_fit_params_img(fh,overlay_mode);

%% contrast panel 
[ph,hglobal_Imin_Info_title,hglobal_Imed_Info_title,hglobal_Imax_Info_title,hglobal_Isd_Info_title,...
    hglobal_Imin_Info_val,hglobal_Imed_Info_val,hglobal_Imax_Info_val,hglobal_Isd_Info_val,...
    hglob_auto,hImintitle,hImin,hImaxtitle,hImax] = ...
    get_contrats_panel_img(fh,strImin,strImax,strImed,strIsd);

%% pointer info 
[hstate,hxtitle,hytitle,hxpos,hypos,...
    hI1title,hI2title,hI1val,hI2val] = get_pointer_info_img(fh,strXpos,strYpos,strI1val,strI2val);


end


function [zoomh,hzoom_title,hzoom_width] = get_zoom_window_img(fh,startzoomwidth)

    set(0,'CurrentFigure',fh);
    zoomh = axes('Parent',fh,'Units','Pixels','Position',[720 25 230 230]);             

    hzoom_title =  uicontrol('Units','pixels',...
                    'Style','text','String','window size',...
                    'Position',[960,150,40,30]);

    hzoom_width =  uicontrol('Units','pixels',...
                    'Style','edit','String',startzoomwidth,...
                    'Position',[960,135,30,15]); 
end

function  [hspots_title,hnumber_of_spots_selected,hnumber_of_spots_selected_title,...
    htotal_number_of_spots,htotal_number_of_spots_title,...
    hspots_xrange_title,hspots_yrange_title,hspots_Irange_title,...
    hspots_xmin,hspots_ymin,hspots_Imin,...
    hspots_xmax,hspots_ymax,hspots_Imax] = get_overlay_spots_settings_img(fh,pos,spots_range,overlay_mode)

set(0,'CurrentFigure',fh);

if ~strcmp(overlay_mode,'spots')
    hspots_title=0; hnumber_of_spots_selected=0; hnumber_of_spots_selected_title=0;
    htotal_number_of_spots=0; htotal_number_of_spots_title=0;
    hspots_xrange_title=0; hspots_yrange_title=0; hspots_Irange_title=0;
    hspots_xmin=0; hspots_ymin=0; hspots_Imin=0;
    hspots_xmax=0; hspots_ymax=0; hspots_Imax=0;
    return
end

hspots_title =  uicontrol('Units','pixels',...
                'Style','text','String','select spots',...
                'Position',[600,230,60,15]);

hnumber_of_spots_selected =  uicontrol('Units','pixels',...
                'Style','text','String',num2str(size(pos,1)),...
                'Position',[570,205,40,15]);            
hnumber_of_spots_selected_title =  uicontrol('Units','pixels',...
                'Style','text','String','spots selected',...
                'Position',[610,205,80,15]);

htotal_number_of_spots =  uicontrol('Units','pixels',...
                'Style','text','String',num2str(size(pos,1)),...
                'Position',[630,190,60,15]);             
htotal_number_of_spots_title =  uicontrol('Units','pixels',...
                'Style','text','String','out of',...
                'Position',[570,190,60,15]);            
            
hspots_xrange_title =  uicontrol('Units','pixels',...
     'Style','text','String','x range',...
     'Position',[600,165,60,15]);
hspots_yrange_title =  uicontrol('Units','pixels',...
    'Style','text','String','y range',...
    'Position',[600,125,60,15]);
hspots_Irange_title =  uicontrol('Units','pixels',...
    'Style','text','String','I range',...
    'Position',[600,45,60,15]);

hspots_xmin =  uicontrol('Units','pixels',...
                'Style','edit','String',num2str(spots_range(1,1)),...
                'Position',[570,150,60,15]);             
hspots_xmax =  uicontrol('Units','pixels',...
                'Style','edit','String',num2str(spots_range(2,1)),...
                'Position',[630,150,60,15]); 
            
hspots_ymin =  uicontrol('Units','pixels',...
                'Style','edit','String',num2str(spots_range(1,2)),...
                'Position',[570,110,60,15]);             
hspots_ymax =  uicontrol('Units','pixels',...
                'Style','edit','String',num2str(spots_range(2,2)),...
                'Position',[630,110,60,15]); 
            
hspots_Imin =  uicontrol('Units','pixels',...
                'Style','edit','String',num2str(spots_range(1,3)),...
                'Position',[570,30,60,15]);             
hspots_Imax =  uicontrol('Units','pixels',...
                'Style','edit','String',num2str(spots_range(2,3)),...
                'Position',[630,30,60,15]); 
            
           
end


function [ph,hglobal_Imin_Info_title,hglobal_Imed_Info_title,hglobal_Imax_Info_title,hglobal_Isd_Info_title,...
    hglobal_Imin_Info_val,hglobal_Imed_Info_val,hglobal_Imax_Info_val,hglobal_Isd_Info_val,...
    hglob_auto,hImintitle,hImin,hImaxtitle,hImax] = ...
    get_contrats_panel_img(fh,strImin,strImax,strImed,strIsd)

    set(0,'CurrentFigure',fh);

    ph = uipanel(fh,'Title','Contrast','Units','pixels',...
                 'Position',[150 600 415 110],'TitlePosition','centertop');

    %info on Image                       
    hglobal_Imin_Info_title = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','text','String','min',...
                    'Position',[0,15,40,15]);            
    hglobal_Imed_Info_title = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','text','String','median',...
                    'Position',[65,15,40,15]);            
    hglobal_Imax_Info_title = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','text','String','max',...
                    'Position',[120,15,40,15]);            
    hglobal_Isd_Info_title = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','text','String','SD',...
                    'Position',[195,15,40,15]);

    hglobal_Imin_Info_val = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','text','String',strImin,...
                    'Position',[00,5,40,15]);            
    hglobal_Imed_Info_val = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','text','String',strImed,...
                    'Position',[65,5,40,15]);            
    hglobal_Imax_Info_val = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','text','String',strImax,...
                    'Position',[120,5,40,15]);                       
    hglobal_Isd_Info_val = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','text','String',strIsd,...
                    'Position',[195,5,40,15]);             

    %contrast settings           
    
    hglob_auto = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','pushbutton',...
                    'String','auto adjust on whole stack',...
                    'Units','pixels','Position',[260,35,150,25]); 

    hImintitle = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','text','String','min Int.',...
                    'Position',[260,17,60,15]);

    hImin = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','edit','String',strImin,...
                    'Position',[260,5,60,15]);         

    hImaxtitle = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','text','String','max Int.',...
                    'Position',[330,17,60,15]);

    hImax = uicontrol('Parent',ph,'Units','pixels',...
                    'Style','edit','String',strImax,...
                    'Position',[330,5,60,15]);                


end

function [hstate,hxtitle,hytitle,hxpos,hypos,...
    hI1title,hI2title,hI1val,hI2val] = get_pointer_info_img(fh,strXpos,strYpos,strI1val,strI2val)

    set(0,'CurrentFigure',fh);
    
    hstate=uicontrol('style','list','String',...
    {'<HTML><FONT COLOR=00FF00>grab</FONT><HTML>',...
    '<HTML><FONT COLOR=FF0000>snap</FONT><HTML>'},...
    'Position',[600,693,60,15]);       

    hxtitle = uicontrol('Style','text','String','x',...
               'Position',[50,645,30,15]);

    hxpos = uicontrol('Style','text','String',strXpos,...
               'Position',[70,645,60,15]);       

    hytitle = uicontrol('Style','text','String','y',...
               'Position',[50,630,60,15]);

    hypos = uicontrol('Style','text','String',strYpos,...
               'Position',[70,630,60,15]);       

    hI1title = uicontrol('Style','text','String','G Int',...
               'Position',[50,615,30,15]);

    hI1val = uicontrol('Style','text','String',strI1val,...
               'Position',[70,615,60,15]);

    hI2title = uicontrol('Style','text','String','R Int 2',...
               'Position',[50,600,30,15]);

    hI2val = uicontrol('Style','text','String',strI2val,...
               'Position',[70,600,60,15]);
end

function [hgauss_title,h_x0_title,h_x0,h_y0_title,h_y0,...
    h_I0_title,h_I0,h_bg0_title,h_bg0,...
    h_sxy0_title,h_sxy0] = get_local_gaussian_fit_params_img(fh,overlay_mode)

    set(0,'CurrentFigure',fh);
    
    if ~strcmp(overlay_mode,'simple')
        hgauss_title=0; h_x0_title=0; h_x0=0; h_y0_title=0; h_y0=0;
        h_I0_title=0; h_I0=0; h_bg0_title=0; h_bg0=0;
        h_sxy0_title=0; h_sxy0=0;
        return;
    end
    
    hgauss_title =  uicontrol('Units','pixels',...
                    'Style','pushbutton','String','local gaussian fit',...
                    'Position',[580,230,100,15]);

    h_x0_title = uicontrol('Units','pixels',...
                    'Style','text','String','x max',...
                    'Position',[580,200,40,15]);
    h_x0 = uicontrol('Units','pixels',...
                    'Style','text','String','undefined',...
                    'Position',[630,200,60,15]);
    h_y0_title = uicontrol('Units','pixels',...
                    'Style','text','String','y max',...
                    'Position',[580,185,40,15]);
    h_y0 = uicontrol('Units','pixels',...
                    'Style','text','String','undefined',...
                    'Position',[630,185,60,15]);            
    h_I0_title = uicontrol('Units','pixels',...
                    'Style','text','String','Int. max',...
                    'Position',[565,150,70,15]);
    h_I0 = uicontrol('Units','pixels',...
                    'Style','text','String','undefined',...
                    'Position',[630,150,60,15]); 
    h_bg0_title = uicontrol('Units','pixels',...
                    'Style','text','String','background',...
                    'Position',[565,135,70,15]);
    h_bg0 = uicontrol('Units','pixels',...
                    'Style','text','String','undefined',...
                    'Position',[630,135,60,15]); 
    h_sxy0_title = uicontrol('Units','pixels',...
                    'Style','text','String','s_xy',...
                    'Position',[580,115,40,15]);
    h_sxy0 = uicontrol('Units','pixels',...
                    'Style','text','String','undefined',...
                    'Position',[630,115,60,15]);     
end


function [xh,yh,hxplot_title,hyplot_title,...
    hxplot_size_txt,hyplot_size_txt,...
    hxplot_range,hyplot_range] = build_axial_view_plots_and_settings_img(fh,xpix,ypix)
    
    set(0,'CurrentFigure',fh);
    xh = axes('Parent',fh,'Units','Pixels','Position',[600 590 330 100]); 
    yh = axes('Parent',fh,'Units','Pixels','Position',[600 440 330 100]); 
    
    hxplot_title =  uicontrol('Units','pixels',...
                    'Style','text','String','profile along x',...
                    'Position',[720,693,100,15]);
    hyplot_title =  uicontrol('Units','pixels',...
                    'Style','text','String','profile along y',...
                    'Position',[720,543,100,15]);
    
    hxplot_size_txt =  uicontrol('Units','pixels',...
                    'Style','text','String','number of voxels to display',...
                    'Position',[940,640,50,45]); 
    hyplot_size_txt =  uicontrol('Units','pixels',...
                    'Style','text','String','number of voxels to display',...
                    'Position',[940,490,50,45]); 
    
    hxplot_range = uicontrol('Units','pixels',...
                    'Style','edit','String',xpix,...
                    'Position',[950,620,30,15]);                  
    hyplot_range = uicontrol('Units','pixels',...
                    'Style','edit','String',ypix,...
                    'Position',[950,470,30,15]);                
end

function [nx,ny,strXpos, strYpos, strI1val, strI2val, ...
    Imin,Imax,Imed,Isd,strImin,strImax,strImed,strIsd,...
    xpix,ypix,startzoomwidth,...
    spots_Imin,spots_Imax,spots_range] = init_variables_for_view_stack_col2_img(stack1,stack2,pos,overlay_mode) 



    %%%%%%%%%%%%%%%%%%%%% variables initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [nx,ny] = size(stack1);
    
    strXpos = '1';     % default x value of the current pixel
    strYpos = '1';     % default y value of the current pixel
    strI1val = num2str(stack1(1,1,1));   % value of the green channel intensity @ the current pixel
    strI2val = num2str(stack2(1,1,1));  % value of the red channel intensity @ the current pixel
    
    %values of global min max and median
    Imin = min(min(stack1));
    Imax = max(max(stack1));
    Imed = median_im(stack1);
    Isd = std_im(stack1);

    strImin = num2str(Imin,'%10.4f');
    strImax = num2str(Imax,'%10.4f');
    strImed = num2str(Imed,'%10.4f');
    strIsd = num2str(Isd,'%10.4f');
    
    %default values of the ranges of the 3 axial views
    xpix = num2str(10);
    ypix = num2str(10);
    
    %default value of the size of the zoom window.
    startzoomwidth = num2str(10);
    if strcmp(overlay_mode,'spots');
        spots_Imin = min(pos(:,3));
        spots_Imax = max(pos(:,3));
        spots_range = [1,1,spots_Imin;nx,ny,spots_Imax];
    else
        spots_Imin = 0; spots_Imax = 0; spots_range = 0;
    end
end