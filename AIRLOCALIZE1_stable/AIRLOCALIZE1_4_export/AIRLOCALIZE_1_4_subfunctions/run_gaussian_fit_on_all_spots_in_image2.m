function detection_result = run_gaussian_fit_on_all_spots_in_image2(varargin)
% version 2: removing display text
tic
verbose = 0;
%% parse arguments
if nargin < 3
    error('not enough arguments to run fit');
else
    rough_pix = varargin{1};
    rough_spots = varargin{2};
    p = varargin{3};
    stack = p.data;
    p = rmfield(p,'data');
    if ndims(stack) == 3
        numdim = 3;
        [nx,ny,nz] = size(stack);
    else
        numdim = 2;
        [nx,ny] = size(stack);
    end
    
    if nargin > 3
        im_num = varargin{4};
        if nargin > 4
            detection_result = varargin{5};
        else
            detection_result = 0;
        end
    else
        im_num = 1;
    end
end

if ~isstruct(detection_result)
    detection_result = struct('rough_spots',{},'final_spots',{},'final_pix',{});
    detection_result(1).rough_spots = cell(1,1);
    detection_result(1).final_spots = cell(1,1);
    detection_result(1).final_pix = cell(1,1);
end
detection_result.rough_spots = rough_spots;


%% initialize arrays
nspots_rough = size(rough_pix,1);       
switch p.fit
    case '3D mask full'
        final_spots = zeros(nspots_rough,6);   % the coordinates output by the fine localization algorithm (in nm)
        final_pix = zeros(nspots_rough,6);   % the coordinates output by the fine localization algorithm (in pix)
        p.cutwidth = [ceil(p.cutsize*p.sigma_xy) ceil(p.cutsize*p.sigma_z)];        
    case '2D mask on local max projection'
        final_spots = zeros(nspots_rough,5);   % the coordinates output by the fine localization algorithm (in nm)
        final_pix = zeros(nspots_rough,5);   % the coordinates output by the fine localization algorithm (in pix)
        p.cutwidth = [ceil(p.cutsize*p.sigma_xy) ceil(p.cutsize*p.sigma_z)];        
    case '3D gaussian fit'
        final_spots = zeros(nspots_rough,8);   % the coordinates output by the fine localization algorithm (in nm)
        final_pix = zeros(nspots_rough,8);   % the coordinates output by the fine localization algorithm (in pix)
        p.cutwidth = [ceil(p.cutsize*p.sigma_xy) ceil(p.cutsize*p.sigma_z)];  
    case '2D gaussian mask'
        final_spots = zeros(nspots_rough,4);   % the coordinates output by the fine localization algorithm (in nm)
        final_pix = zeros(nspots_rough,4);   % the coordinates output by the fine localization algorithm (in pix)
        p.cutwidth = ceil(p.cutsize*p.sigma_xy);          
    case '2D gaussian fit'
        final_spots = zeros(nspots_rough,6);   % the coordinates output by the fine localization algorithm (in nm)
        final_pix = zeros(nspots_rough,6);   % the coordinates output by the fine localization algorithm (in pix)
        p.cutwidth = ceil(p.cutsize*p.sigma_xy); 
end
final_spots(1:nspots_rough,end)=im_num;
final_pix(1:nspots_rough,end)=im_num;
       
%% loop over each pre-detected spot    
for j=1:nspots_rough
    %progress update message
    if(mod(j,100)==0 && verbose==1), disp(['spot ',num2str(j),' out of ' num2str(nspots_rough)]); end

    %Background correction
    switch p.bg_mode
        case 'local plane' 
            [~,stack_bg_corr,new_ctr,ROIlimits] = ...
                gen_linear_interpol_clean_small(...
                stack,rough_pix(j,1:p.numdim),p.cutwidth,p.thickness,'large');
        case 'local median'
            [~,stack_bg_corr,new_ctr,ROIlimits] = ...
                subtract_median(...
                stack,rough_pix(j,1:p.numdim),p.cutwidth,p.thickness,'large','local');
        case 'global median'
            [~,stack_bg_corr,new_ctr,ROIlimits] = ...
                subtract_median(...
                stack,rough_pix(j,1:p.numdim),p.cutwidth,p.thickness,'large','global');
    end
    
    %Gaussian Fit
    switch p.fit
        case '3D mask full'
            % fine localization using the gaussian mask algorithm
            [x0,y0,z0,N0,err0] = gaussian_mask_small(stack_bg_corr,new_ctr,p);
            
            %storing the results
            xmin = ROIlimits(1,1); ymin = ROIlimits(1,2); zmin = ROIlimits(1,3);
            final_pix(j,1)=x0 + xmin - 1;
            final_pix(j,2)=y0 + ymin - 1;
            final_pix(j,3)=z0 + zmin - 1;
            final_pix(j,4)=N0;%Itot
            final_pix(j,5)=err0;%Residuals
            final_spots(j,1:5)=...
                [final_pix(j,1)*p.dx,final_pix(j,2)*p.dy,final_pix(j,3)*p.dz,final_pix(j,4),final_pix(j,5)];
            
        case '2D mask on local max projection'
            %local maxproj
            [img_bg_corr,~,~,~] = local_maxproj(stack_bg_corr,new_ctr,p,'small');

            % fine localization using the gaussian mask algorithm
            [x0,y0,~,N0,err0] = gaussian_mask_small(img_bg_corr,new_ctr(1:2),p);
            
            %storing the results
            xmin = ROIlimits(1,1); ymin = ROIlimits(1,2); zmin = ROIlimits(1,3); 
            final_pix(j,1)=x0 + xmin - 1;
            final_pix(j,2)=y0 + ymin - 1;
            final_pix(j,3)=new_ctr(3) + zmin - 1;
            final_pix(j,4)=N0;  %Itot
            final_pix(j,5)=err0;%Residuals
            final_spots(j,1:5)=...
                [final_pix(j,1)*p.dx,final_pix(j,2)*p.dy,final_pix(j,3)*p.dz,final_pix(j,4),final_pix(j,5)];  
            
        case '3D gaussian fit'
            Gaussout = gaussian_fit_local(stack_bg_corr,new_ctr,p,1);
            xmin = ROIlimits(1,1); ymin = ROIlimits(1,2); zmin = ROIlimits(1,3);
            final_pix(j,1)=Gaussout(5) + xmin - 1;  %xc
            final_pix(j,2)=Gaussout(6) + ymin - 1;  %yc
            final_pix(j,3)=Gaussout(7) + zmin - 1;  %zc
            final_pix(j,4)=Gaussout(8); %Itot
            final_pix(j,5)=Gaussout(9); %Residuals
            
            final_spots(j,1:7)=...
                [final_spots(j,1)*p.dx,final_spots(j,2)*p.dy,final_spots(j,3)*p.dz,final_spots(j,4),...
                final_spots(j,5)*p.dx,final_spots(j,6)*p.dz,final_spots(j,7)];
            
        case '2D gaussian mask'
            % fine localization using the gaussian mask algorithm
            [x0,y0,~,N0,err0] = gaussian_mask_small(stack_bg_corr,new_ctr,p);
            %storing the results
            xmin = ROIlimits(1,1); ymin = ROIlimits(1,2); 
            final_pix(j,1)=x0 + xmin - 1;
            final_pix(j,2)=y0 + ymin - 1;
            final_pix(j,3)=N0;%Itot
            final_pix(j,4)=err0;%Residuals
            final_spots(j,1:3)=...
                [final_pix(j,1)*p.dx,final_pix(j,2)*p.dy,final_pix(j,3)];
            
        case '2D gaussian fit'
            Gaussout = gaussian_fit_local(stack_bg_corr,new_ctr,p,1);
            xmin = ROIlimits(1,1); ymin = ROIlimits(1,2); 
            final_pix(j,1)=Gaussout(4) + xmin - 1;  %xc
            final_pix(j,2)=Gaussout(5) + ymin - 1;  %yc
            final_pix(j,3)=Gaussout(6);  %Itot
            final_pix(j,4)=Gaussout(7);  %residuals
            final_spots(j,1:5)=...
                [final_spots(j,1)*p.dx,final_spots(j,2)*p.dy,final_spots(j,3),...
                final_spots(j,4)*p.dx,final_spots(j,5)];            
    end
  
end 

%% array cleanup & report display
%setting the arrays to the right size, ordering them by decreasing intensity and removing the potential double identifications
if numdim == 3
    [final_spots,final_pix,n_wrong,n_double] = ...
        clean_up_spots_array_clean3(final_pix,[p.dx,p.dy,p.dz],[nx,ny,nz],p.ROIsize,numdim); 
else
    [final_spots,final_pix,n_wrong,n_double] = ...
        clean_up_spots_array_clean3(final_pix,[p.dx,p.dy],[nx,ny],p.ROIsize,numdim); 
end
nspots_found = size(final_spots,1);
detection_result.final_pix = final_pix;
detection_result.final_spots = final_spots;

t = toc;  

disp(['eliminated ' num2str(n_wrong) ' wrong / ' num2str(n_double) ' double identifications; there remain ',num2str(nspots_found),' spots on image ',num2str(im_num)]);
disp(['gaussian mask complete after ', num2str(t), ' s']);
  

end

function [img,zc2,zmin,zmax] = local_maxproj(stack,spot_ctr,p,ROIsize)

zc = spot_ctr(3);
cutwidth_z = p.cutwidth(2);
[~,~,nz]= size(stack);
    
    %the number of planes max projected
    if strcmp(ROIsize,'large')
        zmin = ceil(zc) - floor(cutwidth_z) - floor(p.thickness);
        zmin = max(1,zmin);
        zmax = ceil(zc) + floor(cutwidth_z) + floor(p.thickness);
        zmax = min(nz,zmax);

    else
        zmin = ceil(zc) - floor(cutwidth_z);
        zmin = max(1,zmin);
        zmax = ceil(zc) + floor(cutwidth_z);
        zmax = min(nz,zmax);
    end
    
    zc2 = zc - zmin + 1;
    img = max( stack(:,:,zmin:zmax),[],3 );

end

