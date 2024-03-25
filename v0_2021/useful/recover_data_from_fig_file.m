%retrieves the panels(s) from the current figure
axs = get(gcf, 'Children');

%retrieves the curve(s) from the first panel
pos = get(axs(1), 'Children');

%retrieves the ydata from the first curve
y = get(pos(1),'Ydata');