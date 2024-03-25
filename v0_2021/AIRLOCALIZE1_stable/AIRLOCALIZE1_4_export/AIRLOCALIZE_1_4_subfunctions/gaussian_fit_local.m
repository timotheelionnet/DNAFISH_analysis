function Gaussout = gaussian_fit_local(arr,spot_ctr,p,background_correction)
%% input
%arr is the image or stack
%spots_ctr = [x y z] or [x y] are the guess center coordinates
%p is the structure holding the fit parameters
%p.dx :voxel width
%p.dz : voxel height
%p.sxy : PSF width (same units as that used for the voxel dimensions)
%p.sz : PSF height (same units as that used for the voxel dimensions)
%p.cutsize : range (in PSF width units) of the region around the spot used for fitting. 
%p.type: function type
    %possible values are: 
    %'gaussian' (2D or 3D) - self explanatory; fits background as a constant
    %offset
    %'integrated gaussian' (2D or 3D) - gaussian PSF is integrated over the
    %volume of the pixel/voxel; fits background as a constant
    %offset
    %'integrated gaussian no bg' (2D or 3D) - same as previous except that
    %the background is assumed to be 0
    %'integrated gaussian std z' (for 3D only) - gaussian PSF is integrated over the
    %surface of the xy pixel; std gaussian enveloppe along z; fits background as a constant
    %offset
%p.tol = tolerance of the fit
%p.maxcount = max number of iterations
%background_correction 
    %= 1: subtrat linear interpolation of bakground
    %around spot
    %=0: no background correction
    
   
%% output
%Gaussout is an array with the results:
%(2d) Gaussout= [Imax,bg,sxy,xc,yc,Itot,err]; (no background correction)
%(2d) Gaussout= [Imax,bg,sxy,xc,yc,Itot,err,a,b,c]; (background correction; Ibackground = ax+by+c)
%(3d) Gaussout= [Imax,bg,sxy,sz,xc,yc,zc,Itot,err]; (no background correction)
%(3d) Gaussout= [Imax,bg,sxy,sz,xc,yc,zc,Itot,err,a,b,c,d]; (background correction; Ibackground = ax+by+cz+d)
%Imax is the prefactor of the gaussian (or integrated gaussian)
%bg is the background constant offset
%sxy is the lateral spread of the gaussian
%sz is the axial spread of the gaussian
%xc,yc,zc are the center coordinates
%Itot is the integrated intensity of the spot
%err = sum( (fun(xdata) - ydata)^2 ) / Npoints; where Npoints is the number
%of datapoints

%% retrieving / formatting input data for the fit
numdim = ndims(arr);
xc= spot_ctr(1); 
yc = spot_ctr(2);
if numdim ==3
    zc = spot_ctr(3);
    [nx ny nz] = size(arr);
else
    [nx ny] = size(arr);
end

if ~isstruct(p) %default parameters if
    sxy = double(2.0); 
    sz = double(2.0);
    cutsize = 3;
    tol = 1e-6;
    maxcount = 200;
    
    if numdim == 3
        %p.type = 'integrated gaussian no bg';
        type = 'integrated gaussian std z';
        %p.type = 'standard gaussian';
    else
        type = 'integrated gaussian';
        %p.type = 'standard gaussian';
    end
else
    if isfield(p,'sxy');
        sxy = p.sigma_xy/p.dx;
    else
        sxy = 2.0; 
    end
    
    if isfield(p,'cutsize')
        cutsize = p.cutsize;
    else
        cutsize = 3;
    end
    
    if numdim == 3
        if isfield(p,'sz');
            sz = p.sigma_z/p.dz;
        else
            sz = double(2.0);
        end
    end
    
    if ~isfield(p,'type')
        if numdim == 3
            %p.type = 'integrated gaussian no bg';
            type = 'integrated gaussian std z';
            %p.type = 'standard gaussian';
        else
            type = 'integrated gaussian';
            %p.type = 'standard gaussian';
        end
    else
        if max(strcmp(p.type,{'integrated gaussian no bg';'integrated gaussian std z';...
                'integrated gaussian';'gaussian'})) ~= 0
            type = p.type;
        else
            if numdim == 3
                %p.type = 'integrated gaussian no bg';
                type = 'integrated gaussian std z';
                %p.type = 'standard gaussian';
            else
                type = 'integrated gaussian';
                %p.type = 'standard gaussian';
            end
        end
    end   
    if ~isfield(p,'tol')
        tol = 1e-6;
    else
        tol = p.tol;
    end
    
    if ~isfield(p,'maxcount')
        maxcount = 200;
    else
        maxcount = p.maxcount;
    end
    
    if numdim == 2 && strcmp(p.type,'integrated gaussian std z')
        type = 'integrated gaussian';
    end
end

%% generate ROI around spot to restrict fitting to relevant area
%% and correct for background if necessary
if background_correction == 1
    if numdim == 2
        [~,data,new_ctr,ROIlimits,~,plfit,~] = ...
            gen_linear_interpol_clean_small(arr,[xc,yc],cutsize*sxy,1,'large');
        xmin = ROIlimits(1,1); ymin = ROIlimits(1,2);
        xc2 = new_ctr(1); yc2 = new_ctr(2);

        [nrows, ncols] = size(data);
        [c,r] = meshgrid(1:ncols,1:nrows);
        c = ceil(c) - 0.5;
        r = ceil(r) - 0.5;
        pts = [r c];
        xmax = xmin+size(data,1) - 1;
        ymax = ymin+size(data,2) - 1;
    elseif numdim == 3
        [~,data,new_ctr,ROIlimits,~,plfit,~] = ...
            gen_linear_interpol_clean_small(arr,[xc,yc,zc],[cutsize*sxy,cutsize*sz],1,'large');
        xmin = ROIlimits(1,1); ymin = ROIlimits(1,2); zmin = ROIlimits(1,3);
        xc2 = new_ctr(1); yc2 = new_ctr(2); zc2 = new_ctr(3);
        [nrows, ncols, nups] = size(data);
        [c,r,h] = meshgrid(1:ncols,1:nrows,1:nups);
        c = ceil(c) - 0.5;
        r = ceil(r) - 0.5;
        h = round(h);
        pts = [r c h];
        xmax = xmin+size(data,1) - 1;
        ymax = ymin+size(data,2) - 1;
        zmax = zmin+size(data,3) - 1;
    end
    data = double(data);
else
    if numdim == 3
        npix = [ceil(cutsize*sxy),ceil(cutsize*sxy),ceil(cutsize*sz)];
    else
        npix = [ceil(cutsize*sxy),ceil(cutsize*sxy)];
    end
    xmin = max(ceil(xc)-npix(1),1);
    xmax = min(ceil(xc)+npix(1),nx); 
    ymin = max(ceil(yc)-npix(2),1); 
    ymax = min(ceil(yc)+npix(2),ny);
    if numdim == 3
        zmin = max(ceil(zc)-npix(3),1); 
        zmax = min(ceil(zc)+npix(3),nz); 
        data = arr(xmin:xmax,ymin:ymax,zmin:zmax);
        [nrows, ncols, nups] = size(data);
        [c,r,h] = meshgrid(1:ncols,1:nrows,1:nups);
        c = ceil(c) - 0.5;
        r = ceil(r) - 0.5;
        h = round(h);
        pts = [r c h];
        xc2 = xc - (ceil(xmin) - 0.5);
        yc2 = yc - (ceil(ymin) - 0.5);
        zc2 = zc - round(zmin) +1;
    else
        data = arr(xmin:xmax,ymin:ymax);
        [nrows, ncols] = size(data);
        [c,r] = meshgrid(1:ncols,1:nrows);
        pts = [r c];
        pts = ceil(pts) - 0.5;
        xc2 = xc - (ceil(xmin) - 0.5);
        yc2 = yc - (ceil(ymin) - 0.5);
    end
    data = double(data);
end

%% defining guess values of the fit parameters, their acceptable range, and other fit options
bgmin = double(min(data(:)));
bg = median(double(data(:)));
Imax = double(max(data(:))) - bg;

if numdim == 3
    Coeffs = [  Imax,              bg,            sxy,              sz,             xc2,               yc2,               zc2];
    lb =     [  0,                 bgmin,             0.1,              0.1,            0,                0,                0];
    ub =     [  5*(Imax),       0.5*Imax+bg,      (xmax-xmin)/2,    (zmax-zmin)/2,  (xmax-xmin)+1,      (ymax-ymin)+1,      (zmax-zmin)+1];
else
    Coeffs = [  Imax,              bg,            sxy,              xc2,                yc2];
    lb =     [  0,                 bgmin,             0.1,              0,                 0];
    ub =     [  5*(Imax),       0.5*Imax+bg,      (xmax-xmin)/2,    (xmax-xmin)+1,      (ymax-ymin)+1];
end
options = optimset('TolX',tol,'MaxIter',maxcount,'Display','off');

%% setting the appropriate fitting function
if numdim == 3
    if strcmp(type,'gaussian')
        hfun = @gauss3D;
    elseif strcmp(type,'integrated gaussian')
        hfun = @gauss_integrated3D;
    elseif strcmp(type,'integrated gaussian std z')
        hfun = @gauss_integrated3D_stdZ;
    elseif strcmp(type,'integrated gaussian no bg')
        hfun = @gauss_integrated3D_no_bg;
    end   
elseif numdim ==2
    if strcmp(type,'gaussian')
        hfun = @gauss2D;
    elseif strcmp(type,'integrated gaussian')
        hfun = @gauss_integrated2D;
    elseif strcmp(type,'integrated gaussian no bg')
        hfun = @gauss_integrated2D_no_bg;
    elseif strcmp(type,'integrated gaussian std z')%use integrated PSF if entry was wrong
        hfun = @gauss_integrated2D;
    end    
end

%% actual fit
[Gaussout, resnorm] = lsqcurvefit(hfun, Coeffs, pts,data,lb,ub,options);
if strfind(type,'no bg')
    Gaussout(2) = 0;
end

%% computing the integrated intensity of the fitted spot
if numdim == 3
    switch type
       case 'gaussian'
           xc = Gaussout(5); yc = Gaussout(6); zc = Gaussout(7); sxy = Gaussout(3); sz = Gaussout(4);
           x = (floor(xc - 3*sxy): ceil(xc + 3*sxy )  );
           y = ( floor(yc - 3*sxy): ceil(yc + 3*sxy )  );
           z = ( floor(zc - 3*sz): ceil(zc + 3*sz )  );
           [yy,xx,zz] = meshgrid(y,x,z);
           pts2 = [ceil(xx)-0.5, ceil(yy)-0.5, round(zz)];
           bg = Gaussout(2);
           Gaussout(2) = 0;
           Itot = sum(sum(sum(gauss3D(Gaussout, pts2))));
           Gaussout(2) = bg;
        case 'integrated gaussian'
           I0 = intensity_integrated_gaussian3D(...
                1,1,1,0.5,0.5,1,Gaussout(3),Gaussout(4),1,1,1);
           Itot = 8*Gaussout(1)/I0;
        case 'integrated gaussian std z'
           xc = Gaussout(5); yc = Gaussout(6); zc = Gaussout(7); sxy = Gaussout(3); sz = Gaussout(4);
           x = (floor(xc - 3*sxy): ceil(xc + 3*sxy )  );
           y = ( floor(yc - 3*sxy): ceil(yc + 3*sxy )  );
           z = ( floor(zc - 3*sz): ceil(zc + 3*sz )  );
           [yy,xx,zz] = meshgrid(y,x,z);
           pts2 = [ceil(xx)-0.5, ceil(yy)-0.5, round(zz)];
           bg = Gaussout(2);
           Gaussout(2) = 0;
           Itot = sum(sum(sum(gauss_integrated3D_stdZ(Gaussout, pts2))));
           Gaussout(2) = bg;
        case 'integrated gaussian no bg'
            I0 = intensity_integrated_gaussian3D(...
                1,1,1,0.5,0.5,1,Gaussout(3),Gaussout(4),1,1,1);
            Itot = 8*Gaussout(1)/I0;
    end
    Gaussout(8) = Itot;
    Gaussout(9) = sqrt(resnorm);
    
elseif numdim ==2
    switch type
       case 'gaussian'
           xc = Gaussout(4); yc = Gaussout(5); sxy = Gaussout(3); 
           x = (floor(xc - 3*sxy): ceil(xc + 3*sxy )  );
           y = ( floor(yc - 3*sxy): ceil(yc + 3*sxy )  );
           [yy,xx] = meshgrid(y,x);
           pts2 = [ceil(xx)-0.5, ceil(yy)-0.5];
           bg = Gaussout(2);
           Gaussout(2) = 0;
           Itot = sum(sum(gauss2D(Gaussout, pts2)));
           Gaussout(2) = bg;
        case 'integrated gaussian'
            I0 = intensity_integrated_gaussian2D(...
                1,1,0.5,0.5,Gaussout(3),1,1);
            Itot = 4*Gaussout(1)/I0;
        case 'integrated gaussian no bg'
            I0 = intensity_integrated_gaussian2D(...
                1,1,0.5,0.5,Gaussout(3),1,1);
            Itot = 4*Gaussout(1)/I0;
    end   
    Gaussout(6) = Itot;
    Gaussout(7) = sqrt(resnorm);
end

%% correcting the center position for the the offset of the substack used for
%the fit
if numdim == 3
    Gaussout(5) = Gaussout(5) - 1 + xmin;
    Gaussout(6) = Gaussout(6) - 1 + ymin;
    Gaussout(7) = Gaussout(7) - 1 + zmin;
else
    Gaussout(4) = Gaussout(4) - 1 + xmin;
    Gaussout(5) = Gaussout(5) - 1 + ymin;
end

%% adding the extra parameters from the background equation to the output
if background_correction == 1
    if numdim == 3
        Gaussout(10) = plfit(1);
        Gaussout(11) = plfit(2);
        Gaussout(12) = plfit(3);
        Gaussout(13) = plfit(4);    
    elseif numdim == 2
        Gaussout(8) = plfit(1);
        Gaussout(9) = plfit(2);
        Gaussout(10) = plfit(3);    
    end
end

end


%% the standard gaussians functions formatted for fitting
function I = gauss3D(Coeffs, pts)

%pts = [r c];
%where r is an array of size the image; each value is the row index; c is
%the same for column

%FYI
% Imax = Coeffs(1); 
% bg = Coeffs(2); 
% sxy=Coeffs(3); 
% sz=Coeffs(4); 
% xc=Coeffs(5); 
% yc=Coeffs(6); 
% zc=Coeffs(7);

gridSize = size(pts,2);
r = pts(:,1:gridSize/3,:) ;
c = pts(:,gridSize/3+1:2*gridSize/3,:) ;
h = pts(:,2*gridSize/3+1:gridSize,:) ;

I = Coeffs(2)+ Coeffs(1)*intensity_gaussian3D(r,c,h,...
    Coeffs(5),Coeffs(6),Coeffs(7),Coeffs(3),Coeffs(4),1,1,1);
end

function I = gauss2D(Coeffs, pts)
%pts = [r c];
%where r is an array of size the image; each value is the row index; c is
%the same for column

%FYI
% Imax = Coeffs(1); 
% bg = Coeffs(2); 
% sxy=Coeffs(3); 
% xc=Coeffs(4); 
% yc=Coeffs(5); 

gridSize = size(pts,2);
r = pts( : , 1 : gridSize/2 );
c = pts( : , gridSize/2 + 1 : gridSize );

I = Coeffs(2) + Coeffs(1)*intensity_gaussian2D(r,c,Coeffs(4),Coeffs(5),Coeffs(3),1,1);
end

function I = gauss_integrated3D(Coeffs, pts)
%pts = [r c, h];
%where r is an array of size the stack; each value is the row index; c is
%the same for column, h for planes

%FYI
% Imax = Coeffs(1); 
% bg = Coeffs(2); 
% sxy=Coeffs(3); 
% sz=Coeffs(4); 
% xc=Coeffs(5); 
% yc=Coeffs(6); 
% zc=Coeffs(7);

gridSize = size(pts,2);
r = pts(:,1:gridSize/3,:);
c = pts(:,gridSize/3+1:2*gridSize/3,:);
h = pts(:,2*gridSize/3+1:gridSize,:);

%compute max intensity  so prefactor (Coeffs(1)) is actually
%Imax
I0 = intensity_integrated_gaussian3D(...
     1,1,1,0.5,0.5,1,Coeffs(3),Coeffs(4),1,1,1);

%compute intensity over all datapoints
I = Coeffs(2) + Coeffs(1)*intensity_integrated_gaussian3D...
    (r,c,h,Coeffs(5),Coeffs(6),Coeffs(7),Coeffs(3),Coeffs(4),1,1,1)/I0;

end

function I = gauss_integrated3D_stdZ(Coeffs, pts)
%pts = [r c, h];
%where r is an array of size the stack; each value is the row index; c is
%the same for column, h for planes

%FYI
% Imax = Coeffs(1); 
% bg = Coeffs(2); 
% sxy=Coeffs(3); 
% sz=Coeffs(4); 
% xc=Coeffs(5); 
% yc=Coeffs(6); 
% zc=Coeffs(7);

gridSize = size(pts,2);
r = pts(:,1:gridSize/3,:);
c = pts(:,gridSize/3+1:2*gridSize/3,:);
h = pts(:,2*gridSize/3+1:gridSize,:);

%compute max intensity so prefactor (Coeffs(1)) is actually
%Imax
I0 = intensity_integrated_gaussian3D_stdZ(...
    1,1,1,0.5,0.5,1,Coeffs(3),Coeffs(4),1,1,1);

%compute intensity over all datapoints
I = Coeffs(2) + Coeffs(1)*intensity_integrated_gaussian3D_stdZ...
    (r,c,h,Coeffs(5),Coeffs(6),Coeffs(7),Coeffs(3),Coeffs(4),1,1,1)/I0;

end

function I = gauss_integrated3D_no_bg(Coeffs, pts)
%pts = [r c, h];
%where r is an array of size the stack; each value is the row index; c is
%the same for column, h for planes

%FYI
% Imax = Coeffs(1); 
% bg = Coeffs(2) (useless in this function)
% sxy=Coeffs(3); 
% sz=Coeffs(4); 
% xc=Coeffs(5); 
% yc=Coeffs(6); 
% zc=Coeffs(7);

gridSize = size(pts,2);
r = pts(:,1:gridSize/3,:);
c = pts(:,gridSize/3+1:2*gridSize/3,:);
h = pts(:,2*gridSize/3+1:gridSize,:);

%compute max intensity  so prefactor (Coeffs(1)) is actually
%Imax
I0 = intensity_integrated_gaussian3D(...
     1,1,1,0.5,0.5,1,Coeffs(3),Coeffs(4),1,1,1);

%compute intensity over all datapoints
I = Coeffs(1)*intensity_integrated_gaussian3D...
    (r,c,h,Coeffs(5),Coeffs(6),Coeffs(7),Coeffs(3),Coeffs(4),1,1,1)/I0;

end

function I = gauss_integrated2D_no_bg(Coeffs, pts)
%pts = [r c];
%where r is an array of size the image; each value is the row index; c is
%the same for column

%FYI
% Imax = Coeffs(1); 
% bg = Coeffs(2) (useless in this function)
% sxy=Coeffs(3); 
% xc=Coeffs(4); 
% yc=Coeffs(5); 

gridSize = size(pts,2);
r = pts( : , 1 : gridSize/2 );
c = pts( : , gridSize/2 + 1 : gridSize );

%compute max intensity  so prefactor (Coeffs(1)) is actually
%Imax
I0 = intensity_integrated_gaussian2D(...
     1,1,0.5,0.5,Coeffs(3),1,1);

%compute intensity over all datapoints
I = Coeffs(1)*intensity_integrated_gaussian2D...
    (r,c,Coeffs(4),Coeffs(5),Coeffs(3),1,1)/I0;

end

function I = gauss_integrated2D(Coeffs, pts)
%pts = [r c];
%where r is an array of size the image; each value is the row index; c is
%the same for column

%FYI
% Imax = Coeffs(1); 
% bg = Coeffs(2) 
% sxy=Coeffs(3); 
% xc=Coeffs(4); 
% yc=Coeffs(5); 

gridSize = size(pts,2);
r = pts( : , 1 : gridSize/2 );
c = pts( : , gridSize/2 + 1 : gridSize );

%compute max intensity  so prefactor (Coeffs(1)) is actually
%Imax
I0 = intensity_integrated_gaussian2D(...
     1,1,0.5,0.5,Coeffs(3),1,1);


%compute intensity over all datapoints
I = Coeffs(2) + Coeffs(1)*intensity_integrated_gaussian2D...
    (r,c,Coeffs(4),Coeffs(5),Coeffs(3),1,1)/I0;


end
