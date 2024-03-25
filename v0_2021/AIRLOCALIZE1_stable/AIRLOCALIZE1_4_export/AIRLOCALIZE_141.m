function varargout = AIRLOCALIZE_141(varargin)
%no argument: you load a tiff image yourself or point the program to the directory where your images are stored.

%one (optional) argument: it should be either: a 2D image or a 3D stack 
%second (optional) argument is a structure p that contains the parameters for
%detection (no user intervention needed)

persistent params
tStart = tic;
%params: structure composed by the input parameters for detection:

AIRLOCALIZE_version = 1;

notice_str = sprintf(['Algorithm For Three Dimensional Localization of Single mRNAs and Transcription Sites in Cells and Tissues\n',...
'© 2010 Albert Einstein College of Medicine of Yeshiva University\n',... 
'All rights reserved.\n\n',...
'Image Loader\n',...
'Copyright (c) 2008, Andriy Nych\n',...
'All rights reserved.\n\n',...
'Tiff reader\n',...
'Copyright (c) 2006, Francois Nedelec\n',...
'All rights reserved. \n\n',...
'Version 1.4\n\n',...
'Please reference Lionnet et al, Nature Methods 2011\n',...
'Questions:\n',...
'Timothee Lionnet, lionnett@janelia.hhmi.org']);
 
license_window = dispwin('copyright notice', notice_str);

 
%% checking the type of argument / parameter initialization

%this function sets all parameters for the detection program
params = set_detection_parameters_141(params,AIRLOCALIZE_version,varargin);
if strcmp(params.mode,'cancel'), clear('mode','tStart'); return; end

%% detection
if strcmp(params.mode, 'single file input')
    
    %retrieve raw data and attach it to params structure
    params = retrieve_raw_data(params);
    fh = dispwin('Progress','Processing Image...');

    if params.dense_spots == 0
        [detection_result,params] = ...
            perform_detection_on_single_image_once(params);
    elseif params.dense_spots ==1
        [detection_result,params] = ...
            perform_detection_on_single_image_once_dense_spots(params);
    end
    
    %cleanup
    if isfield(params,'data'), params = rmfield(params,'data'); end
    if isfield(params,'smooth'), params = rmfield(params,'smooth'); end
    
    detection_result.params = params;    
    varargout{1} = detection_result;
    
    %% saving the results and parameters
    if ~isempty(params.save_dirname)
        [~, fname_no_ext] = fileparts(params.data_stackname);
        if params.numdim == 2
            loc3 = [detection_result.final_pix(:,1:3),zeros(size(detection_result.final_pix,1),1)];
        elseif params.numdim == 3
            loc3 = [detection_result.final_pix(:,1:4),zeros(size(detection_result.final_pix,1),1)];
        end

        %ascii list of the spots (coordinates in pix)
        save(fullfile(params.save_dirname,[fname_no_ext,'.loc3']),'loc3','-ascii');
        %ascii list of the parameters
        save_structure_as_text(params,fullfile(params.save_dirname,[fname_no_ext,'.par3']));

        %binary matlab file holding both the results and parameters
        save(fullfile(params.save_dirname,[fname_no_ext,'.det']),'detection_result');

        %% generating the image that contains the spots
        if isfield(params,'output_spot_image')
            if params.output_spot_image == 1
                fh = dispwin('Progress','generating output image...',fh);

                if params.AIRLOCALIZE_version == 1
                    %output an image / stack where the center of the spots are
                    %marked by squares
                    if params.numdim == 3
                            out1 = convert_coordinates_to_stack3(detection_result.final_pix(:,1:3),params.nx,params.ny,params.nz,'ValMode','ones','Size',params.sigma_xy);
                    elseif params.numdim == 2
                            out1 = convert_coordinates_to_stack3(detection_result.final_pix(:,1:2),params.nx,params.ny,'ValMode','ones','Size',params.sigma_xy);
                    end

                else
                    %output an image where the objects footprints are marked by an
                    %area which value is equal to the number of mRNAs

                    out1 = generate_image_of_objects(detection_result.cc,params);
                end
                [~, fname_no_ext] = fileparts(params.data_stackname);
                save_as_tiff(uint16(out1),fullfile(params.save_dirname,[fname_no_ext,'_spots.tif']));
                clear('out1');

            end
        end    
    end
    close(fh);
elseif strcmp(params.mode, 'directory input')
    [detection_result,params]  = perform_detection_from_dir(params);
    varargout{1} = detection_result;
    
elseif strcmp(params.mode, 'batch movies input')
    [detection_result,params]  = perform_detection_batch_movies_from_dir(params);
    varargout{1} = detection_result;
       
elseif  strcmp(params.mode,'movie directory input') 
    [detection_result,params] = perform_detection_movie_from_dir(params);
    varargout{1} = detection_result;
    
elseif strcmp(params.mode,'movie single file input')
    [detection_result,params] = perform_detection_movie_from_file(params);
    varargout{1} = detection_result;
    
elseif (strcmp(mode,'failed'))
    disp('data input failed');
    varargout{1} = 0;
    return 
else
    disp('incorrect file input mode');
    varargout{1} = 0;
    return
end

t = toc(tStart);

disp(['detection over after ' num2str(t) ' s']);
clear('stack','mode','detection_result_movie','t','tStart');
close(license_window);

end

function [detection_result,p] = perform_detection_on_single_image_once_dense_spots(varargin)
%takes 2d or 3d images as input (called stack in both cases)
%perform_detection_on_single_image_once(stack,p)
%perform_detection_on_single_image_once(stack,p,'image number', i,'smoothed stack', smooth)

%% parsing input
%parameters
p = varargin{1};

%img number
im_num=1; %default value for the image number
for nn = 1:nargin
    if strcmp(varargin{nn},'image number')
        if nargin>nn
            im_num = varargin{nn+1};
        end
    end
end

%loading the stacks if it is not already stored in the parameters structure
p = retrieve_raw_data(p);

% retrieve the detected objects (or detect the current image)
p = retrieve_detected_objects(p);

stack = p.data;
p = rmfield(p,'data');

%% detection/quantification
if strcmp(p.dense_mode,'exclusive')
   [detection_result,p] = detect_gaussian_mask_exclusive3(stack,p,im_num); 
elseif strcmp(p.dense_mode,'inclusive')
   [detection_result,p] = detect_gaussian_mask_inclusive4(stack,p,im_num);
end

clear('stack');

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

function [detection_result,p] = perform_detection_on_single_image_once(varargin)
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

function [detection_result,p] = perform_detection_batch_movies_from_dir(p)

if ~isfield(p,'exclusion_string')
    p.exclusion_string = {[]};
end
if ~isfield(p,'inclusion_string')
    p.inclusion_string = {[]};
end

flist = get_clean_file_list(p.data_dirname,p.inclusion_string,...
    p.exclusion_string,strcmp(p.dirmode,'recursive'),0);
flist = flist(cellfun(@is_recognized_img_format,flist));

if strcmp(p.save_dirname,p.data_dirname)
    save_in_same_dir = 1;
else
    save_in_same_dir = 0;
end

for i=1:numel(flist)
        if save_in_same_dir ==1
            [savedir{i} , ~, ~ ] = fileparts(flist{i});
        else
            savedir{i} = p.save_dirname;
        end
end

i = 1; %file index in the directory
j = 1; %number of recognized image files
fh = dispwin('Progress','Starting Batch Process...');
while i <= numel(flist)
    
    if is_recognized_img_format(flist{i})
        if j == 1
            p.data_stackname = cell(1);
        end
        %p.data_stackname{j,1} = fullfile(list(i).dir,list(i).name);
        p.data_stackname{j,1} = flist{i};
        p.curfname = p.data_stackname{j,1};
        
        fh = dispwin('Progress',['opening file ',num2str(i),'/',num2str(numel(flist)),': ',p.data_stackname{j,1}],fh);
        
        p.data = double(timtiffread( p.data_stackname{j,1} ));
        tmpname = p.data_stackname;
        p.data_stackname = p.data_stackname{j,1};
        [detection_result{j},p] = perform_detection_movie_from_file(p);   
        p.data_stackname = tmpname;
        j = j+1;
    end
    i = i +1;
end
clear('i','j','list');


if isfield(p,'global_nx')
    p = rmfield(p,'global_nx');
end
if isfield(p,'global_ny')
    p = rmfield(p,'global_ny');
end
if p.numdim == 3
    if isfield(p,'global_nz')
        p = rmfield(p,'global_nz');
    end
end

% house keeping
clear('smooth','stack','list','i','j','rough_spots','rough_pix','filename');
close(fh);
end


function [detection_result,p] = perform_detection_from_dir(p)

if ~isfield(p,'exclusion_string')
    p.exclusion_string = {[]};
end
if ~isfield(p,'inclusion_string')
    p.inclusion_string = {[]};
end

flist = get_clean_file_list(p.data_dirname,p.inclusion_string, p.exclusion_string,strcmp(p.dirmode,'recursive'),0);
flist = flist(cellfun(@is_recognized_img_format,flist));

if strcmp(p.save_dirname,p.data_dirname)
    save_in_same_dir = 1;
else
    save_in_same_dir = 0;
end

for i=1:numel(flist)
        if save_in_same_dir ==1
            [savedir{i} , ~, ~ ] = fileparts(flist{i});
        else
            savedir{i} = p.save_dirname;
        end
end

i = 1; %file index in the directory
j = 1; %number of recognized image files
fh = dispwin('Progress','Starting Batch Process...');
while i <= numel(flist)
    
    if is_recognized_img_format(flist{i})
        if j == 1
            p.data_stackname = cell(1);
        end
        %p.data_stackname{j,1} = fullfile(list(i).dir,list(i).name);
        p.data_stackname{j,1} = flist{i};
        p.curfname = p.data_stackname{j,1};
        
        fh = dispwin('Progress',['opening file ',num2str(i),'/',num2str(numel(flist)),': ',p.data_stackname{j,1}],fh);
        if p.dense_spots == 0
            [detection_result{j},p] = perform_detection_on_single_image_once(p);
        elseif p.dense_spots ==1    
            [detection_result{j},p] = perform_detection_on_single_image_once_dense_spots(p);
        end

        detection_result{j}.params = p;
        [~, fname_no_ext] = fileparts(p.data_stackname{j,1});
        
        %ascii list of the spots (coordinates in pix)
        if p.numdim == 2
            loc3 = [detection_result{j}.final_pix(:,1:3),zeros(size(detection_result{j}.final_pix,1),1)];
        elseif p.numdim == 3
            loc3 = [detection_result{j}.final_pix(:,1:4),zeros(size(detection_result{j}.final_pix,1),1)];
        end
        save(fullfile(savedir{i},[fname_no_ext,'.loc3']),'loc3','-ascii');
        
        %ascii list of the parameters
        tmpname = p.data_stackname;
        p.data_stackname = p.data_stackname{j,1};
        
        %saving the size of the full image
        if p.numdim ==3
            if j == 1
                p.global_nx = cell(1);
                p.global_ny = cell(1);
                p.global_nz = cell(1);
            end
            [p.global_nx{j,1}, p.global_ny{j,1}, p.global_nz{j,1}]= deal(p.nx,p.ny,p.nz);
        elseif p.numdim ==2
            if j == 1
                p.global_nx = cell(1);
                p.global_ny = cell(1);
            end
            [p.global_nx{j,1}, p.global_ny{j,1}]= deal(p.nx,p.ny);
        end
        
        if p.numdim == 2
            [global_nx,global_ny] = deal(p.global_nx,p.global_ny);
            p = rmfield(p,'global_nx');
            p = rmfield(p,'global_ny');
        elseif p.numdim == 3
            [global_nx,global_ny,global_nz] = deal(p.global_nx,p.global_ny,p.global_nz);
            p = rmfield(p,'global_nx');
            p = rmfield(p,'global_ny');
            p = rmfield(p,'global_nz');
        end
        if isfield(p,'inclusion_string')
            p = rmfield(p,'inclusion_string');
        end
        if isfield(p,'exclusion_string')
            p = rmfield(p,'exclusion_string');
        end
        save_structure_as_text(p,fullfile(savedir{i},[fname_no_ext,'.par3']));
        p.data_stackname = tmpname;
        
        
        if p.numdim == 2
            [p.global_nx,p.global_ny] = deal(global_nx,global_ny);
        elseif p.numdim == 3
            [p.global_nx,p.global_ny,p.global_nz] = deal(global_nx,global_ny,global_nz);
        end
        %binary matlab file holding the results (and parameters)
        save(fullfile(savedir{i},[fname_no_ext,'.det']),'detection_result');
        
        
        %% generating the image that contains the spots
        if isfield(p,'output_spot_image')
            if p.output_spot_image == 1

                disp('generating output image...\n');
                if p.numdim == 3
                        out1 = convert_coordinates_to_stack3(detection_result{j}.final_pix(:,1:3),p.nx,p.ny,p.nz,'ValMode','ones','Size',p.sigma_xy);
                elseif p.numdim == 2
                        out1 = convert_coordinates_to_stack3(detection_result{j}.final_pix(:,1:2),p.nx,p.ny,'ValMode','ones','Size',p.sigma_xy);
                end
                save_as_tiff(uint16(out1),fullfile(savedir{i},[fname_no_ext,'_spots.tif']));
                clear('out1');
            end
        end

        j = j+1;
    end
    i = i +1;
end
clear('i','j','list');


if isfield(p,'global_nx')
    p = rmfield(p,'global_nx');
end
if isfield(p,'global_ny')
    p = rmfield(p,'global_ny');
end
if p.numdim == 3
    if isfield(p,'global_nz')
        p = rmfield(p,'global_nz');
    end
end

% house keeping
clear('smooth','stack','list','i','j','rough_spots','rough_pix','filename');
close(fh);
end

function [detection_result,p] = perform_detection_movie_from_dir(p)

if ~isfield(p,'exclusion_string')
    p.exclusion_string = {[]};
end
if ~isfield(p,'inclusion_string')
    p.inclusion_string = {[]};
end

flist = get_clean_file_list(p.data_dirname,p.inclusion_string, p.exclusion_string,strcmp(p.dirmode,'recursive'),0);
flist = flist(cellfun(@is_recognized_img_format,flist));  
fh = dispwin('Progress','Starting Movie Analysis...');
i = 1; %file index in the directory
j = 1; %number of recognized image files

while i <= size(flist,1)
    
    if is_recognized_img_format(flist{i})
        
        if j == 1
            p.data_stackname = cell(1);
        end
        p.data_stackname{j,1} = flist{i};
        p.curfname = p.data_stackname{j,1};
        
        stack = timtiffread(p.data_stackname{j,1});
        
        %saving the size of the image
        if p.numdim == 3
            if j == 1
                p.nx = cell(1);
                p.ny = cell(1);
                p.nz = cell(1);
            end
            [p.nx{j,1}, p.ny{j,1}, p.nz{j,1}]= size(stack);
        elseif p.numdim ==2
            if j == 1
                p.nx = cell(1);
                p.ny = cell(1);
            end
            [p.nx{j,1}, p.ny{j,1}]= size(stack);
        end
        
        %cropping the image if ROI has been selected
        if p.select_ROI == 1
            if p.numdim == 3
                if ~(size(stack,1) == p.Xmax - p.Xmin+1 && size(stack,2) == p.Ymax - p.Ymin+1 && size(stack,3) == p.Zmax - p.Zmin+1)
                   stack = stack(p.Xmin:p.Xmax,p.Ymin:p.Ymax,p.Zmin:p.Zmax);
                end
            elseif p.numdim == 2
                if ~(size(stack,1) == p.Xmax - p.Xmin+1 && size(stack,2) == p.Ymax - p.Ymin+1)
                   stack = stack(p.Xmin:p.Xmax,p.Ymin:p.Ymax);
                end
            end
        end
        p.data = stack;
        clear('stack');
        
        fh = dispwin('Progress',['opening file ',num2str(i),'/',num2str(numel(flist)),': ',p.data_stackname{j,1}],fh);
        if p.dense_spots == 0
            [detection_result{j},p] = perform_detection_on_single_image_once(p,'image number', j);
        elseif p.dense_spots ==1    
            [detection_result{j},p] = perform_detection_on_single_image_once_dense_spots(p,'image number', j);
        end
        
        if j== 1
            if p.numdim == 2
                loc3 = [detection_result{j}.final_pix(:,1:3),(j-1)*ones(size(detection_result{j}.final_pix,1),1)];
            elseif p.numdim == 3
                loc3 = [detection_result{j}.final_pix(:,1:4),(j-1)*ones(size(detection_result{j}.final_pix,1),1)];
            end
        else
            if p.numdim == 2
                loc3 = [loc3;detection_result{j}.final_pix(:,1:3),(j-1)*ones(size(detection_result{j}.final_pix,1),1)];
            elseif p.numdim == 3
                loc3 = [loc3;detection_result{j}.final_pix(:,1:4),(j-1)*ones(size(detection_result{j}.final_pix,1),1)];
            end
        end
        
        %% generating the image that contains the spots
        if isfield(p,'output_spot_image')
            if p.output_spot_image == 1

                disp('generating output image...\n');
                if p.numdim == 3
                        out1 = convert_coordinates_to_stack3(detection_result{j}.final_pix(:,1:3),p.nx{j,1},p.ny{j,1},p.nz{j,1},'ValMode','ones','Size',p.sigma_xy);
                elseif p.numdim == 2
                        out1 = convert_coordinates_to_stack3(detection_result{j}.final_pix(:,1:2),p.nx{j,1},p.ny{j,1},'ValMode','ones','Size',p.sigma_xy);
                end
                [~, fname_no_ext] = fileparts(p.data_stackname{j,1});
                save_as_tiff(uint16(out1),fullfile(p.save_dirname,[fname_no_ext,'_spots.tif']));
                clear('out1');
            end
        end
        
        j = j+1;
    end
    i = i +1;
end
clear('i','j','flist');

% output formatted for tracking routine called utrack 
detection_result{1}.movieInfo = output_trajectories_as_structure_from_res2(detection_result);
if isfield(p,'data')
    p = rmfield(p,'data');
end
if isfield(p,'smooth')
    p = rmfield(p,'smooth');
end
detection_result{1}.params = p;
save(fullfile(p.save_dirname,'batch_spot_detection.loc3'),'loc3','-ascii');             
save(fullfile(p.save_dirname,'batch_spot_detection.det'),'detection_result');
p.nx = p.nx{1,1};
p.ny = p.ny{1,1};
if p.numdim == 3
    p.nz = p.nz{1,1};
end
        
%saving the size of the full image

if isfield(p,'inclusion_string')
    p = rmfield(p,'inclusion_string');
end
if isfield(p,'exclusion_string')
    p = rmfield(p,'exclusion_string');
end
if isfield(p,'data_stackname')
    %p.flist_str = sprintf('%s\\n',p.data_stackname{:});
    %p = rmfield(p,'data_stackname');
end
            
save_structure_as_text(p,fullfile(p.save_dirname,'batch_spot_detection.par3'));

% house keeping
clear('smooth','stack','list','i','j','rough_spots','rough_pix','filename');
close(fh);
end

function [detection_result,p] = perform_detection_movie_from_file(p)
 
%single file containing a 2D movie (i.e. 3D stack)
if isfield(p,'data')
    if ndims(p.data) == 3
        stack = p.data;
        p = rmfield(p,'data');
    else
        stack = timtiffread(fullfile( p.data_dirname , p.data_stackname ));
        p = rmfield(p,'data');
    end
else
    stack = timtiffread(fullfile( p.data_dirname , p.data_stackname ));
end 

%saving the size of the image
[p.nx, p.ny,~]= size(stack);

%cropping the movie is ROI has been selected
if p.select_ROI == 1
    if ~(size(stack,1) == p.Xmax - p.Xmin+1 && size(stack,2) == p.Ymax - p.Ymin+1 && size(stack,3) == p.Zmax - p.Zmin+1)
         stack = stack(p.Xmin:p.Xmax,p.Ymin:p.Ymax,p.Zmin:p.Zmax);
    end
end 
fh = dispwin('Progress','Starting Movie Analysis...');
for j = 1:size(stack,3)
    
    fh = dispwin('Progress',['Processing frame ',num2str(j),'/',num2str(size(stack,3))],fh);
    p.data = stack(:,:,j);
    [detection_result{j},p] = perform_detection_on_single_image_once(...
        p,'image number', j);

    if j== 1
        if p.numdim == 2
            loc3 = [detection_result{j}.final_pix(:,1:3),(j-1)*ones(size(detection_result{j}.final_pix,1),1)];
        elseif p.numdim == 3
            loc3 = [detection_result{j}.final_pix(:,1:4),(j-1)*ones(size(detection_result{j}.final_pix,1),1)];
        end
    else
        if p.numdim == 2
            loc3 = [loc3;detection_result{j}.final_pix(:,1:3),(j-1)*ones(size(detection_result{j}.final_pix,1),1)];
        elseif p.numdim == 3
            loc3 = [loc3;detection_result{j}.final_pix(:,1:4),(j-1)*ones(size(detection_result{j}.final_pix,1),1)];
        end
    end
    
    %% generating the image that contains the spots
    if isfield(p,'output_spot_image')
        if p.output_spot_image == 1

            disp('generating output image...\n');
            out1 = convert_coordinates_to_stack3(detection_result{j}.final_pix(:,1:2),p.nx,p.ny,'ValMode','ones','Size',p.sigma_xy);
            [~, fname_no_ext] = fileparts(p.data_stackname);
            imwrite(uint16(out1),fullfile(p.save_dirname,[fname_no_ext,'_spots.tif']),'Compression','none','WriteMode','append');
            clear('out1');
        end
    end
end

clear('i','j','list');

% output formatted for tracking routine called utrack
detection_result{1}.movieInfo = output_trajectories_as_structure_from_res2(detection_result);
    
%binary matlab file holding both the results and parameters
if isfield(p,'data')
    p = rmfield(p,'data');
end
if isfield(p,'smooth')
    p = rmfield(p,'smooth');
end

detection_result{1}.params = p;
[~, fname_no_ext] = fileparts(p.data_stackname);
save(fullfile(p.save_dirname,[fname_no_ext,'.det']),'detection_result');
save_structure_as_text(p,fullfile(p.save_dirname,[fname_no_ext,'.par3']));
save(fullfile(p.save_dirname,[fname_no_ext,'.loc3']),'loc3','-ascii');         

% house keeping
clear('smooth','stack','list','i','j','rough_spots','rough_pix','filename');
close(fh);
end


