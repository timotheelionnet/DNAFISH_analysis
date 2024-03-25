function [ xpos,xvals ] = find_tick_values_for_img_axes(xmin,xmax,npts,zerocentered,Nopt)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

xpos = [];
xvals = [];
%find ticks for image data with values

if zerocentered
    if abs(xmax) ~= abs(xmin)
        disp('the data is not symmetric, cannot center on zero');
        zerocentered = 0;
    end
    
    if mod(npts,2) == 0
        disp('even number of pixels, cannot center on zero');
        zerocentered = 0;
    end
end



    
if zerocentered    
    nhalf = (npts-1)/2; 
    Nopt = Nopt/2;
    
    nsteps = 1:2*Nopt;
    pixperstep = (1:nhalf)';
    
    x = repmat(nsteps,numel(pixperstep),1) .* repmat(pixperstep,1,numel(nsteps));
    
    x(x > nhalf) = -Inf;
    
    [x,ix] = max(x);
    
    ix = ix(~isinf(x));
    x = x(~isinf(x));
    
    [~,jxopt] = min( abs(x - nhalf)/nhalf + abs(ix-Nopt)/Nopt );
    
    nsteps = ix(jxopt);
    pixperstep = jxopt;
    nmax = x(jxopt);
    
    xpos =  (nhalf - nmax) +1 : pixperstep : 2*nsteps*pixperstep + (nhalf - nmax) +1;
    
    dx = (xmax-xmin)/(npts-1);
    
    xvals = (-nsteps*pixperstep : pixperstep : nsteps*pixperstep) * dx;
else
    
    nsteps = 1:2*Nopt;
    pixperstep = (1:npts)';
    
    x = repmat(nsteps,numel(pixperstep),1) .* repmat(pixperstep,1,numel(nsteps)) +1;
    
    x(x > npts) = -Inf;
    
    [x,ix] = max(x);
    
    ix = ix(~isinf(x));
    x = x(~isinf(x));
    
    [~,jxopt] = min( abs(x - npts)/npts + abs(ix-Nopt)/Nopt );
    
    nsteps = ix(jxopt);
    pixperstep = jxopt;
    nmax = x(jxopt);
    
    xpos = 1 : pixperstep : nsteps*pixperstep +1;
    
    dx = (xmax-xmin)/(npts-1);
    
    xvals = (xpos-1) * dx + xmin;
    
    
end





end

