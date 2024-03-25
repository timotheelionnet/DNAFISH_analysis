function [xc,yc] = compute_radial_hist_circular_density( x,y )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

binsize = x(2) - x(1);

xc = x;

areanorm = pi*(xc+binsize/2).^2 - pi*(max(0,xc-binsize/2)).^2;

yc = y./areanorm;






end

