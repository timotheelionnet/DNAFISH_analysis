function params = retrieve_raw_data(params)
% function that retrieves the raw data image/stack

%% Case 1: data is already loaded in the parameters structure
if isfield(params,'data')
    if ndims(params.data) == params.numdim && max(size(params.data)) > 1  
        params = crop_image(params);
        return
    end
end
    
%% Case 2: data needs to be loaded
%current file name is defined: load file directly
if isfield(params,'curfname')
    if ~strcmp(params.curfname,'null')
        params.data = timtiffread(params.curfname);
    else
        params.data = timtiffread(fullfile(params.data_dirname,params.data_stackname));   
    end
    
    if params.numdim == 3
        [params.nx,params.ny,params.nz] = size(params.data);
    elseif params.numdim == 2
        [params.nx,params.ny] = size(params.data);
    end
    params = crop_image(params);
    return
end

%current file name is undefined: load open file GUI
[fname,sourcedir,fidx] = uigetfile('*.tif;*.stk;*.lsm','Pick Image File Used to Select Parameters',...
    params.data_dirname);
if fidx == 0, params.mode = 'cancel'; return; end

params.curfname = fullfile(sourcedir,fname);
params.data = timtiffread(params.curfname);

if strcmp(params.mode,'movie single file input') ...
        strcmp(params.mode,'batch movies input')
    params.data = params.data(:,:,params.movie_frame_used_to_select_parameters);
end

if params.numdim == 3
    [params.nx,params.ny,params.nz] = size(params.data);
elseif params.numdim == 2
    [params.nx,params.ny] = size(params.data);
end
params = crop_image(params);

end

function params = crop_image(params)
%croping the image if an ROI was selected

if params.select_ROI == 1
    if isfield(params,'Xmax')
       if ~(size(params.data,1) == params.Xmax - params.Xmin+1 && size(params.data,2) == params.Ymax - params.Ymin+1) 
           if params.numdim == 3
               if ~(size(params.data,3) == params.Zmax - params.Zmin+1)
                    params.data = params.data(params.Xmin:params.Xmax,params.Ymin:params.Ymax,params.Zmin:params.Zmax);
               end
           elseif params.numdim ==2
               params.data = params.data(params.Xmin:params.Xmax,params.Ymin:params.Ymax);
           end
       end
    end    
end
    
end

