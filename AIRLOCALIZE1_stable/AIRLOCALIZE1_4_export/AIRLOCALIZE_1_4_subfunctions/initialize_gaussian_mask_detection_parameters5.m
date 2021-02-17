function p = initialize_gaussian_mask_detection_parameters5(p)
    
%% structure that collects all the parameters of the particle detection

if ~isstruct(p)
    p = struct(...
    'sigma_xy',{},'sigma_z',{},...
    'dx',{},'dy',{},'dz',{},...
    'cutsize',{},'tol',{},'thickness',{},...
    'thresh',{},...
    'max_spots',{},'ROIsize',{},'method',{},'fit',{},...
    'maxcount',{},...
    'data_fname',{},'data_stackname',{},'data_dirname',{},...
    'save_dirname',{});
    p(1).numdim = 3;
    
    % see below for a description of the parameters. these are some
    % defaults for a given microscope
    p.sigma_xy = 2;
    p.sigma_z = 2;
    p.dx = 64;
    p.dy = 64;
    p.dz = 200;
    p.cutsize = 3;
    p.tol = 0.01;
    p.thickness = 1;
    p.thresh(1).level = 3;
    p.thresh.units = 'SD';
    p.thresh.sd = 0;
    p.filter.nlo = 2.15;
    p.filter.nhi = 1;
    %p.filter.numdim = 2;
    %p.filter.width = 0.01; 
    p.max_spots = 200000; 
    p.ROIsize = 2;
    p.fit = '3d mask';
    p.maxcount = 100;
    p.data_fname = [];
    p.data_stackname = [];
    p.data_dirname = [];
    p.save_dirname = [];
    p.type = 'integrated gaussian';
    p.set_psf_manually = 1;
    p.set_thresh_manually = 1;
    p.data = 0;
    p.smooth = 0;
    p.bg_mode = 'local plane';
    p.mode = 'null';
end

if ~isfield(p,'numdim')
    p(1).numdim = 3;
end
%--------------------------------------------------------------------------

%% these are the default parameters
    
%PSF width in pixels in resp. (xy) and z
if ~isfield(p,'sigma_xy')
    p.sigma_xy = 2; 
end
if p.numdim == 3
    if ~isfield(p,'sigma_z')
        p.sigma_z = 2;
    end
end

%voxel size in nm
if ~isfield(p,'dx')
    p.dx = 64;
end

if ~isfield(p,'dy')
    p.dy = 64;
end

if p.numdim == 3
    if ~isfield(p,'dz')
        p.dz = 200; 
    end
end
if ~isfield(p,'cutsize')
    p.cutsize = 3;      %size of the ROI surrounding each local intensity maximum 
end                    %over which the gaussian mask is performed (in PSF width units)

if ~isfield(p,'tol')
p.tol = 0.01;           %tolerance for the convergence of the gaussian mask
end

if ~isfield(p,'thickness')
    p.thickness = 1;    %size of the region surrounding the 3D-rectangle (see above) 
end                    %that the program uses to compute the local background

if ~isfield(p,'thresh')
    p.thresh(1).fc_xy = 1.5;    %spatial wavelength (in PSF units)
    if p.numdim == 3             %used for low cutoff in the bandpass filter of the image               
        p.thresh.fc_z = 1.5;        %spatial wavelength (in PSF units)
                                %used for low cutoff in the bandpass filter of the image
    end                            
    p.thresh.level = 3;         %value used to threshold the pixels on the filtered image 
                                %(either in absolute or SD units)
    p.thresh.units = 'SD';      %= 'SD' : 'level' parameter is in units of the SD of the filtered image
                                %= 'absolute' : 'level' parameter is in
                                %absolute units 
    p.thresh.sd = 0;            %value of the SD of the filtered image
end                           


if ~isfield(p,'max_spots')
p.max_spots = 200000;         %maximum number of detected spots on the image 
end
if ~isfield(p,'ROIsize')
    p.ROIsize = 2;              %minimum distance allowed between two distinct spots
end

if ~isfield(p,'threshold')
    p.method = 'threshold';     %method of predection = 'threshold' or ='spottiness' (obsolete)
end

if ~isfield(p,'fit')
    p.fit = '3d mask';
    %p.fit = '2d mask';         %local maxproj
    %p.fit = '3d fit';          %local 3d gaussian NLLS fit
end

if ~isfield(p,'maxcount')
    p.maxcount = 100;           %max number of iterations of the gaussian mask
end
if ~isfield(p,'data_fname')
    p.data_fname = [];          %name of the original data file (if applicable) 
end
if ~isfield(p,'data_stackname')
    p.data_stackname = [];      %name of the variable holding the data
end
if ~isfield(p,'data_dirname')
    p.data_dirname = [];        %name of the directory from which the data files are extracted (if applicable)
end
if ~isfield(p,'save_dirname')
    p.save_fname = 'D:/temp/';          %name of the file in which this structure is stored
end
if ~isfield(p,'type')
    p.type = 'integrated gaussian';
end
if ~isfield(p,'set_psf_manually')
    p.set_psf_manually = 1;
end
if ~isfield(p,'set_thresh_manually')
    p.set_thresh_manually = 1;
end
if ~isfield(p,'data')
    p.data = 0;
end
if ~isfield(p,'smooth')
    p.smooth = 0;
end
 %-------------------------------------------------------------------------   
end