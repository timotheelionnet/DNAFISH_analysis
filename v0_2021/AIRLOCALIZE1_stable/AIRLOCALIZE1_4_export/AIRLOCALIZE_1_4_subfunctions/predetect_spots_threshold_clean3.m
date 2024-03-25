function [rough_spots,rough_pix] = predetect_spots_threshold_clean3(p,i)
%takes an image (smoothimage) and given a set of parameters (hold in
%structure p), returns the positions and intensities of pixels that are predetected based on a threshold. i
%is the image number (used for metadata purpose only)
%spots positions are given with the 'subpixel' convention i.e. center of
%the bottom left pixel = (0.5,0.5,0.5)

%threshold parameters
%p.thresh.level: the value of the threshold (might get affected by a multiplying factor, see below)
%p.thresh.units: can be either 'absolute' or 'SD'. If 'absolute' option, the
%image is threhsolded using the value p.thresh.level. If 'SD' option, the
%threshold value is equal to p.thresh.level*p.thresh.sd

%once the pixels/voxels above threhsold have been selected, a filter is run
%so that if neighbouring pixels/voxels have values above the threshold,
%only the pixel/voxel w/ the maximum intensity within p.ROIsize pixels.


%% threshold level computation
%getting the threshold level that will be used to select spots based on the
%pixel intensities on the filtered image
if strcmp(p.thresh.units,'absolute')
    threshInt = p.thresh.level;
    disp(['threshold value is ' num2str(threshInt) ' in absolute units']);
elseif strcmp(p.thresh.units,'SD')
    if p.numdim == 3
        p.thresh.sd = std_stack(p.smooth);
    elseif p.numdim ==2
        p.thresh.sd = std_im(p.smooth);
    end
    threshInt = p.thresh.level*p.thresh.sd;
    disp(['threshold value is ' num2str(p.thresh.level) ' in SD units, i.e. ' num2str(threshInt) ' in absolute units']);
elseif strcmp(p.thresh.units,'automatic')
    threshInt = p.thresh.thr_min + p.thresh.level*(p.thresh.thr_opt - p.thresh.thr_min);
    disp(['threshold value is ' num2str(p.thresh.level) ' in auto units, i.e. ' num2str(threshInt) ' in absolute units']);
    
end




%% finding the points above the threshold in the filtered image 
%then keep only the local maxima
maxima = find_isolated_maxima_clean2(p.smooth,threshInt,1);
p = rmfield(p,'smooth');

if isempty(maxima)
    msg = 'predetected no spots;';
    disp(msg);
    rough_spots = [];
    rough_pix = [];
    return
end

%ordering the maxima by descending intensity value
if p.numdim == 3
    maxima = sort_array_by_col_value(maxima,4,'descend');
elseif p.numdim == 2
    maxima = sort_array_by_col_value(maxima,3,'descend');
end

%truncating the array if it has more points than allowed
if size(maxima,1) > p.max_spots
    maxima = maxima(1:p.max_spots,:);
    msg = ['predetected ',num2str(size(maxima,1)),' spots, reducing their number to the max allowed number = ', num2str(p.max_spots)];
    disp(msg);
else
    msg = ['predetected ',num2str(size(maxima,1)),' spots;'];
    disp(msg);
end

%--------------------------------------------------------------------------

%% building the arrays holding the lists of maxima coordinates and values
nspots = size(maxima,1);
if p.numdim == 3
    rough_spots = zeros(nspots,5); 
    rough_spots(1:nspots,1) = maxima(1:nspots,1)*p.dx - 0.5*p.dx;    
    rough_spots(1:nspots,2) = maxima(1:nspots,2)*p.dy - 0.5*p.dy;    
    rough_spots(1:nspots,3) = maxima(1:nspots,3)*p.dz - 0.5*p.dz; 
    rough_spots(1:nspots,4) = maxima(1:nspots,4);
    rough_spots(1:nspots,5) = i;

    rough_pix(1:nspots,1:4) = maxima(1:nspots,1:4);
    rough_pix(1:nspots,5) = i;
elseif p.numdim == 2
    rough_spots = zeros(nspots,4); 
    rough_spots(1:nspots,1) = maxima(1:nspots,1)*p.dx - 0.5*p.dx;    
    rough_spots(1:nspots,2) = maxima(1:nspots,2)*p.dy - 0.5*p.dy;    
    rough_spots(1:nspots,3) = maxima(1:nspots,3);
    rough_spots(1:nspots,4) = i;

    rough_pix(1:nspots,1:3) = maxima(1:nspots,1:3);
    rough_pix(1:nspots,4) = i;
end
%--------------------------------------------------------------------------

%% housekeeping
clear('maxima','msg','n_spots','threshInt');

end


function [sort_arr,indices] = sort_array_by_col_value(arr,ncol,orientation)
%arr is the array, 
%ncol is the number of the column which values decide of the ordering of
%the rows
%orientation is a string, either 'ascend' or 'descend'

ny = size(arr,2);
if ny<ncol
    disp('your array does not have enough columns');
end

[~, indices] = sort(arr(:,ncol),orientation);
nx = size(arr,1);
sort_arr(1:nx,1:ny) = arr(indices(1:nx),1:ny);
clear('nx','ny','sortVal');
end

