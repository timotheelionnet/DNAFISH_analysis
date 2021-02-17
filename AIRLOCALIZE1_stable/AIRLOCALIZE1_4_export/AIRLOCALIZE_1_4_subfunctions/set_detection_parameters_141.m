function params = set_detection_parameters_141(params,AIRLOCALIZE_version,varargin)


if ~isstruct(params) 
    params = initialize_gaussian_mask_detection_parameters5(params);
else
    params = reset_parameters(params);
end
params.AIRLOCALIZE_version = AIRLOCALIZE_version;

if nargin >=1
    need_data_type_input =0;
    optarg = varargin{1};
    if ~isempty(optarg)
        if ischar(optarg{1})
            if exist(optarg{1},'file')
                params.data = timtiffread(optarg{1});
                params.curfname = optarg{1};
                [params.data_dirname,params.data_stackname,ext] = fileparts(optarg{1});
                params.data_stackname = [params.data_stackname,ext];
                params.save_dirname = params.data_dirname; 
                if ndims(params.data) == 2
                    params.numdim = 2;
                    params.mode = 'single file input';
                    [params.nx,params.ny] = size(params.data);
                elseif ndims(params.data) == 3
                    params.numdim = 3;
                    params.mode = 'single file input';
                    [params.nx,params.ny,params.nz] = size(params.data);
                else
                    need_data_type_input =1;
                end 
            else
                need_data_type_input =1;
            end
        elseif isnumeric(optarg{1})
            params.curfname = 'variable';
            params.data_dirname = [];
            params.data_stackname = 'variable';
            params.save_dirname = []; 
            if ndims(optarg{1}) == 2
                params.numdim = 2;
                params.mode = 'single file input';
                params.data = optarg{1};
                params.thresh.sd = std(params.data(:));
            elseif ndims(optarg{1}) == 3
                params.numdim = 3;
                params.mode = 'single file input';
                params.data = optarg{1};
                params.thresh.sd = std(params.data(:));
            else
                need_data_type_input =1;
            end 
        else
            need_data_type_input =1;
        end
    else
        need_data_type_input =1;
    end
else
    need_data_type_input =1;
end

    

%% choose the mode (single file vs. directory)
if need_data_type_input
    params = set_detection_mode_141(params);
    if strcmp(params.mode,'cancel'), return; end
end

%% select image file / directory
if need_data_type_input
    params = get_image_file_location(params);
end
if strcmp(params.mode,'movie single file input') ...
        || strcmp(params.mode,'batch movies input')
    params.movie_frame_used_to_select_parameters = 1; 
end
if strcmp(params.mode,'cancel'), return; end

%% set parameters
params = detection_parameters_interface5(params);
if strcmp(params.mode,'cancel'), return; end

if params.dense_spots == 1 && params.set_thresh_manually == 0
    params = set_dense_spots_parameters(params);
end
if strcmp(params.mode,'cancel'), return; end

%% select ROI (optional)
if isfield(params,'select_ROI')
    if params.select_ROI == 1
        params = select_ROI(params);
    end
end
if strcmp(params.mode,'cancel'), return; end

%% set PSF width manually (optional)
if params.set_psf_manually == 1
    params = set_psf_size_manually(params);
end
if strcmp(params.mode,'cancel'), return; end

%% set Int threshold manually (optional) 
%with or without Morphology filters

%automatic threshold optimization
if strcmp(params.thresh.units,'automatic')
    disp('computing auto threshold...');
    params = predict_threshold(params);   
end

if params.set_thresh_manually == 1 
   params = set_threshold_manually(params); 
end

if strcmp(params.mode,'cancel'), return; end

%in the case of batch processing, I clear the data attached to the
%params structure
if ~strcmp(params.mode,'single file input')
    if isfield(params,'data'), params = rmfield(params,'data'); end
    if isfield(params,'smooth'), params = rmfield(params,'smooth'); end
    if isfield(params,'curfname'), params = rmfield(params,'curfname');  end  
end


end


%% subfunctions
function params = set_morphology_manually(params)

    params = retrieve_smoothed_data(params);
    params = explore_stack_2channels_morphology2(params); 

end

function params = get_image_file_location(params)
if strcmp(params.mode,'single file input')
    [fname,sourcedir,fidx] = uigetfile('*.tif;*.stk;*.lsm','Select Source Image File');
    if fidx == 0, params.mode = 'cancel'; return; end
    params.data_dirname = sourcedir;
    params.data_stackname = fname;
    params.save_dirname = params.data_dirname;
    params.curfname = fullfile(sourcedir,fname);
    
elseif strcmp(params.mode,'directory input')
    sourcedir = uigetdir('','Select Source Images Directory');
    if sourcedir == 0, params.mode = 'cancel'; return; end
    params.save_dirname = sourcedir;
    params.data_dirname = sourcedir;
    
elseif strcmp(params.mode,'movie directory input')
    sourcedir = uigetdir('','Select Source Images Directory');
    if sourcedir == 0, params.mode = 'cancel'; return; end
    params.save_dirname = sourcedir;
    params.data_dirname = sourcedir;
elseif strcmp(params.mode,'batch movies input')
    sourcedir = uigetdir('','Select Source Images Directory');
    if sourcedir == 0, params.mode = 'cancel'; return; end
    params.save_dirname = sourcedir;
    params.data_dirname = sourcedir;

elseif strcmp(params.mode,'movie single file input')
    [fname,sourcedir,fidx] = uigetfile('*.tif;*.stk;*.lsm','Select Source Image File');
    if fidx == 0, params.mode = 'cancel'; return; end
    params.data_dirname = sourcedir;
    params.data_stackname = fname;
    params.save_dirname = params.data_dirname;
    params.curfname = fullfile(sourcedir,fname);
    
end

clear('fname','sourcedir','fidx');
end

function params = set_psf_size_manually(params)

    params = retrieve_raw_data(params);
    
    %reducing to a 2D frame in caase of a movie
    if strcmp(params.mode,'movie single file input') ...
            || strcmp(params.mode,'batch movies input')
        params.data = params.data(:,:,params.movie_frame_used_to_select_parameters);
    end
    
    if params.numdim ==3
        p = rmfield(params, 'data');
        res = explore_stack_v2(params.data,p);
    else
        p = rmfield(params, 'data');
        res = explore_img_v2(params.data,p);
    end
    
    if max(res(:)) == 0
        params.mode = 'cancel';
        disp('you did not choose any spot. ABORT!');
        return
    end
    
    if params.numdim == 3
        params.sigma_xy = mean(res(:,6));
        params.sigma_z = mean(res(:,7));
    elseif params.numdim == 2    
        params.sigma_xy = mean(res(:,5));
    end

end

function params = set_threshold_manually(params)

    params = retrieve_smoothed_data(params);
    maxsmooth = max(params.smooth(:))
    
    if params.numdim ==3
        params = explore_stack_2channels2(params); 
    elseif params.numdim ==2
        params = explore_img_2channels2_41(params); 
    end
    
end

function params = select_ROI(params)
    
params = retrieve_raw_data(params);
if params.numdim ==2 && ~strcmp(params.mode,'movie single file input') ...
        && ~strcmp(params.mode,'batch movies input')
    res = TimSelectROI(params.data);
else
    res = TimSelectROI3D(params.data);
end

if isstruct(res)
    %inverted coordinate convention
    params.Xmin = res.Ymin;
    params.Ymin = res.Xmin;

    params.Xmax = res.Ymax;
    params.Ymax = res.Xmax;

    if params.numdim == 3 || strcmp(params.mode,'movie single file input')...
            || strcmp(params.mode,'batch movies input')
        params.Zmin = res.Zmin;
        params.Zmax = res.Zmax;
    end
else
    params.select_ROI = 0;
end

if isfield(params,'data')
    if ndims(params.data) == 3
        params.data = params.data(params.Xmin:params.Xmax,params.Ymin:params.Ymax,params.Zmin:params.Zmax);
    else
        params.data = params.data(params.Xmin:params.Xmax,params.Ymin:params.Ymax);
    end
end

clear('img','stack','res');

end

