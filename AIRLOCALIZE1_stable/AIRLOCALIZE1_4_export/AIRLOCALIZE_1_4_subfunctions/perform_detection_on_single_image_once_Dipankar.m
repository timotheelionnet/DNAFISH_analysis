function [detection_result,p] = perform_detection_on_single_image_once_Dipankar(varargin)
%takes 2d or 3d images as input (called stack in both cases)
%perform_detection_on_single_image_once(p)
%perform_detection_on_single_image_once(p,'image number', i,'smoothed stack', smooth)

%% parsing input
im_num=1; %default value for the image number
p = varargin{1};

for nn = 1:nargin
    if strcmp(varargin{nn},'image number')
        if nargin>nn
            im_num = varargin{nn+1};
        end
    end
end

%retrieve smoothed data
p = retrieve_smoothed_data(p);

%% predetection on one image
[rough_spots,rough_pix] = predetect_spots_threshold_clean3(p,im_num);    

%% detection/quantification

detection_result = run_gaussian_fit_on_all_spots_in_image(rough_pix,rough_spots,p,im_num);
detection_result = rmfield(detection_result,'rough_spots');

%% correcting for the position offset in the case of an ROI
if isfield(p,'select_ROI')
      if p.select_ROI == 1
          detection_result.final_pix(:,1) = detection_result.final_pix(:,1) + p.Xmin-1;
          detection_result.final_spots(:,1) = detection_result.final_spots(:,1) + (p.Xmin-1)*p.dx;
          detection_result.final_pix(:,2) = detection_result.final_pix(:,2) + p.Ymin-1;
          detection_result.final_spots(:,2) = detection_result.final_spots(:,2) + (p.Ymin-1)*p.dx;
      end
end

%% housekeeping
if isfield(p,'data'), p = rmfield(p,'data'); end
if isfield(p,'smooth'), p = rmfield(p,'smooth'); end
clear('rough_pix','rough_spots','im_num');

end