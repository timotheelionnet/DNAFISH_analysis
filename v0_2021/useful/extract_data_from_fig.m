function [xdata,ydata] = extract_data_from_fig(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% h=gcf;
% axesObjs = get(h, 'Children');  %axes handles
axesObjs = gca;
dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
objTypes = get(dataObjs, 'Type');  %type of low-level graphics object

xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs, 'YData');

end

