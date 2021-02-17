function [fh,ah] = create_fig_with_axis_controls(title_msg,varargin)

%[fh,ah] = create_fig_with_axis_controls(title_msg)
%[fh,ah] = create_fig_with_axis_controls(title_msg,limits)
%[fh,ah] = create_fig_with_axis_controls(title_msg,X,Y)
%[fh,ah] = create_fig_with_axis_controls(title_msg,X,Y,limits)
%[fh,ah] = create_fig_with_axis_controls(title_msg,X,Y,limits,color)
%color is either a string or a cell array (number of elements identical as
%number of X columns; each element is eith a linespec argument (e.g. '-or')
%or a 1x3 RGB set of values between 0 and 1.

%X and Y should be 2D arrays of identical size. The program plots colums
%separately as individual curves.

%collecting optional arguments
if ~isempty(varargin)
    %limits(1) = xmin 
    %limits(2) = xmax 
    %limits(3) = ymin 
    %limits(4) = ymax
    if numel(varargin) == 1
        limits = varargin{1};
        xdata = [];
        ydata = [];
        color = set_auto_color(size(xdata,2));
        
    elseif numel(varargin) == 2
        xdata = varargin{1};
        ydata = varargin{2};
        if ndims(xdata)~=2 || ndims(ydata)~=2
            fh = 0;
            ah = 0;
            dipswin('Plot Error','X and Y data should be 2D arrays');
            return
        end
        if sum(size(xdata)== size(ydata))~=2
            fh = 0;
            ah = 0;
            dipswin('Plot Error','X and Y data have different sizes');
            return
        end
        limits = [min(xdata(:)), max(xdata(:)),min(ydata(:)), max(ydata(:))];
        color = set_auto_color(size(xdata,2));
        
    elseif numel(varargin) == 3
        xdata = varargin{1};
        ydata = varargin{2};
        limits = varargin{3};
        color = set_auto_color(size(xdata,2));
        
    elseif numel(varargin) == 4
        xdata = varargin{1};
        ydata = varargin{2};
        limits = varargin{3};
        color = varargin{4};
        if iscell(color)
            if numel(color) ~= size(xdata,2)
                dipswin('Plot Warning','Could Not Read Color Command; Using Auto Colors');
                color = set_auto_color(size(xdata,2));
            end
        else
            if ischar(color)
                color = {color};
                color = repmat(color,1,size(xdata,2));
            elseif isdouble(color)
                color = {color};
                color = repmat(color,1,size(xdata,2));
            else
                dipswin('Plot Warning','Could Not Read Color Command; Using Auto Colors');
                color = set_auto_color(size(xdata,2));
            end
        end
    end    
else
    xdata = [];
    ydata = [];
    limits = [0,1,0,1];
end

%setting the figure up
fh = figure('Units','pixels',...
                  'NumberTitle','off',...
                  'Name',title_msg,...
                  'Position',[200 200 600 400],...
                  'Visible','off'); 

%main plot
ah = axes('Parent',fh,'Units','normalize','Position',[0.2,0.2,0.6,0.6],'Tag','ah');
%plot data
if ~ishold(ah)
    hold(ah);
end
if ~isempty(xdata)
    for i=1:size(xdata,2)
        if ischar(color{i})
            plot(xdata(:,i),ydata(:,i),color{i},'DisplayName',['data column ',num2str(i)]);
        else
            plot(xdata(:,i),ydata(:,i),'DisplayName',['data column ',num2str(i)],'Color',color{i});
        end
    end
end

%enforce limits
xlim([limits(1),limits(2)]);
ylim([limits(3),limits(4)]);

%build axis controls
set_axis_limits_controls(limits,fh);

%Callbacks
h = guihandles(fh);
set(h.hxmin,'Callback', {@change_axes_limits,fh,ah});
set(h.hxmax,'Callback', {@change_axes_limits,fh,ah});
set(h.hymin,'Callback', {@change_axes_limits,fh,ah});
set(h.hymax,'Callback', {@change_axes_limits,fh,ah});
set(h.hauto,'Callback', {@change_axes_limits,fh,ah});


set(fh,'Visible','on');
drawnow;
end

function change_axes_limits(src,eventdata,fh,ah)
    srctag= get(src,'Tag');
    limits = [get(ah,'Xlim'),get(ah,'Ylim')];
    h = guihandles(fh);
    switch srctag
        case 'hxmin'
            limits(1) = str2double(get(src,'String'));
            if isnan(limits(1))
                limits(1) = limits(2) - 1;
                set(src,'String',num2str(limits(1)));
            elseif limits(1) >= limits(2)
                limits(1) = limits(2) - 1;
                set(src,'String',num2str(limits(1)));
            end
        case 'hxmax'
            limits(2) = str2double(get(src,'String'));
            if isnan(limits(2))
                limits(2) = limits(1) + 1;
                set(src,'String',num2str(limits(2)));
            elseif limits(1) >= limits(2)
                limits(2) = limits(1) + 1;
                set(src,'String',num2str(limits(2)));
            end
        case 'hymin'
            limits(3) = str2double(get(src,'String'));
            if isnan(limits(3))
                limits(3) = limits(4) - 1;
                set(src,'String',num2str(limits(3)));
            elseif limits(3) >= limits(4)
                limits(3) = limits(4) - 1;
                set(src,'String',num2str(limits(3)));
            end
        case 'hymax'
            limits(4) = str2double(get(src,'String'));
            if isnan(limits(4))
                limits(4) = limits(3) + 1;
                set(src,'String',num2str(limits(4)));
            elseif limits(3) >= limits(4)
                limits(4) = limits(3) + 1;
                set(src,'String',num2str(limits(4)));
            end
        case 'hauto'
            lh=findall(ah,'type','line');
            xdata = get(lh,'xdata');
            ydata = get(lh,'ydata');                
            for i=1:numel(lh)
                curxdata = xdata{i};
                curydata = ydata{i};
                if i==1
                    limits = [min(curxdata(:)), max(curxdata(:)), min(curydata(:)), max(curydata(:))];
                else
                    limits(1) = min([limits(1),min(curxdata(:))]);
                    limits(2) = max([limits(2),max(curxdata(:))]);
                    limits(3) = min([limits(3),min(curydata(:))]);
                    limits(4) = max([limits(4),max(curydata(:))]);
                end    
            end 
            set(h.hxmin,'String',num2str(limits(1)));
            set(h.hxmax,'String',num2str(limits(2)));
            set(h.hymin,'String',num2str(limits(3)));
            set(h.hymax,'String',num2str(limits(4)));
    end
    set(ah,'XLim',limits(1:2));
    set(ah,'YLim',limits(3:4));
    drawnow;
end

function set_axis_limits_controls(limits,fh)

%Xmin label
htxt = uicontrol('Parent',fh,...
    'Units','normalized',...
    'Style','text','String','Xmin: ',...
    'Position',[0.2 0.05 0.1 0.1]);

set(htxt,'units','character');              
pos = get(htxt,'Position');
set(htxt,'Position',[pos(1), pos(2) ,6,1.5]);
%Xmin Control
hedit = uicontrol('Parent',fh,...
    'Units','normalized',...
    'Tag','hxmin',...
    'Style','edit','String',num2str(limits(1)),...
    'Position',[0.2 0.05 0.1 0.1]);

set(hedit,'units','character');              
set(hedit,'Position',[pos(1)+6, pos(2) ,6,1.5]);


%Xmax label
htxt = uicontrol('Parent',fh,...
    'Units','normalized',...
    'Style','text','String','Xmax: ',...
    'Position',[0.8 0.05 0.1 0.1]);

set(htxt,'units','character');              
pos = get(htxt,'Position');
set(htxt,'Position',[pos(1)-12, pos(2) ,6,1.5]);
%Xmax Control
hedit = uicontrol('Parent',fh,...
    'Units','normalized',...
    'Tag','hxmax',...
    'Style','edit','String',num2str(limits(2)),...
    'Position',[0.2 0.05 0.1 0.1]);

set(hedit,'units','character');              
set(hedit,'Position',[pos(1)-6, pos(2) ,6,1.5]);


%Ymin label
htxt = uicontrol('Parent',fh,...
    'Units','normalized',...
    'Style','text','String','Ymin: ',...
    'Position',[0.05 0.2 0.1 0.1]);

set(htxt,'units','character');              
pos = get(htxt,'Position');
set(htxt,'Position',[pos(1), pos(2) ,6,1.5]);
%Ymin Control
hedit = uicontrol('Parent',fh,...
    'Units','normalized',...
    'Tag','hymin',...
    'Style','edit','String',num2str(limits(3)),...
    'Position',[0.1 0.2 0.1 0.1]);

set(hedit,'units','character');              
set(hedit,'Position',[pos(1)+6, pos(2) ,6,1.5]);


%Ymax label
htxt = uicontrol('Parent',fh,...
    'Units','normalized',...
    'Style','text','String','Ymax: ',...
    'Position',[0.05 0.8 0.1 0.1]);

set(htxt,'units','character');              
pos = get(htxt,'Position');
set(htxt,'Position',[pos(1), pos(2)-1.5 ,7,1.5]);
%Xmin Control
hedit = uicontrol('Parent',fh,...
    'Units','normalized',...
    'Tag','hymax',...
    'Style','edit','String',num2str(limits(4)),...
    'Position',[0.2 0.1 0.1 0.1]);

set(hedit,'units','character');              
set(hedit,'Position',[pos(1)+7, pos(2)-1.5 ,6,1.5]);


%Auto Limits
hauto = uicontrol('Parent',fh,...
    'Units','normalized',...
    'Tag','hauto',...
    'Style','pushbutton',...
    'String','auto',...
    'Position',[0.45 0.05 0.1 0.1]);

set(hauto,'units','character'); 
pos = get(hauto,'Position');
set(hauto,'Position',[pos(1), pos(2) ,6,1.5]);

end

function color = set_auto_color(n)

    RGB = [1,0,0;...
     0,0.5,1;...
     0.1,0.8,0;...
     1,0.9,0;...
     0,0,0;...
     
     1,0.8,0.8;...
     0.5,0.85,1;...
     0.6,0.2,0;...
     0.8,1,0.6;...
     0.5,0.5,0.5;...
     
     1,0.6,0.2;...
     0,0.9,0.9;...
     0.1,1,0.6;...
     1,0,1;...
     0.2,0.2,0.2;...
     ];

%     idx = 0:0.2:1;
%     npts = numel(idx);
%     R = idx(mod(1:n,npts)+1);
%     G = idx(mod(ceil((1:n)/2),npts)+1);
%     B = idx(mod(ceil((1:n)/3),npts)+1);
    
    npts = size(RGB,1);
    color = RGB(rem(1:n,npts)+1,:);
    color = mat2cell(color,ones(1,n),3);
end



