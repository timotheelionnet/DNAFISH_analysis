function params = retrieve_smoothed_data2(params)
% function that retrieves the raw data image/stack
% version 2: calls smooth_image_and_subtract_background6 (vs 5) where there
% is no display window 
%% Case 1: smoothed data is already loaded in the parameters structure
if isfield(params,'smooth')
    if ndims(params.smooth) == params.numdim && size(params.smooth,params.numdim) ~= 1
        return
    end
end
    

%% Case 2: smoothed data needs to be loaded

%retrieve raw data
params = retrieve_raw_data(params);

%reducing to a 2D frame in case of a movie
if strcmp(params.mode,'movie single file input')
    params.data = params.data(:,:,params.movie_frame_used_to_select_parameters);
end

%smooth raw data
params = smooth_image_and_subtract_background6(params);

end