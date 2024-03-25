%compute_Intensity_Spot_Size_vs_dist

%this is the address of the loc file I generated with Analyze_1color_FISH_distribution_for_fig.m
locfile_signal_only = '/Users/lionnett/Documents/Salina Long/two_color/101716/C1-comp_tim2c_2to1_r1s6_rh_all_recon_trim_signal_spots.loc3';
%briefly, I
%1) ran AIRLOCALIZE on the original data stack (parameters in script above)
tifdata = '/Users/lionnett/Documents/Salina Long/two_color/101716/C1-comp_tim2c_2to1_r1s6_rh_all_recon_trim.tif';

%2) converted pixels in nm (88 nm/pix in x-y, 100 nm/pix in z)
%3) flipped x and y (for ViSP - Careful here I need to reflip!)
%4) added a dummy frame column (for ViSP)

%voxel dimensions in nm
xvox = 88;
zvox = 100;

%parameter structure (used by AIRLOCALIZE subfunctions)
p.cutsize = 3;
p.sigma_xy = 0.8;
p.sigma_z = 1.5;
p.cutwidth = [ceil(p.cutsize*p.sigma_xy) ceil(p.cutsize*p.sigma_z)]; 
p.thickness = 1;


%% load spot localization results
rnm = dlmread(locfile_signal_only);
rnm(:,1:2) = [rnm(:,2),rnm(:,1)]; %re-flip the x and y!
rpix = rnm;
rpix(:,1:2) = rpix(:,1:2)/xvox;
rpix(:,3) = rpix(:,3)/zvox;

%% load image dataset
stack = timtiffread(tifdata);

%% loop through spots and perform a gaussian fit with free sigmas as well as intensities
rgauss = zeros(size(rpix,1),6);

%spots that I found to generate non-linear least square fitting issues
%(triggering an error) and I just ignore here
%problem_idx = [1553,2418,3496,4277,4686,5102,5528,5714,6039,6041,6100,6101];

for i=1:size(rpix,1);                                           
    if mod(i,100)==0
        disp(['spot ',num2str(i),'/',num2str(size(rpix,1))]);
    end
    
    if ~ismember(i,problem_idx)
    %generate background corrected spot
    [~,stack_bg_corr,new_ctr,ROIlimits] = ...
                gen_linear_interpol_clean_small(...
                stack,rpix(i,1:3),p.cutwidth,p.thickness,'large');
    
    %fit data to a 3D gaussian
    Gaussout = gaussian_fit_local(stack_bg_corr,new_ctr,p,0);
            xmin = ROIlimits(1,1); 
            ymin = ROIlimits(1,2); 
            zmin = ROIlimits(1,3);
            rgauss(i,1)=Gaussout(5) + xmin - 1;  %xc
            rgauss(i,2)=Gaussout(6) + ymin - 1;  %yc
            rgauss(i,3)=Gaussout(7) + zmin - 1;  %zc
            rgauss(i,4)=Gaussout(8); %Itot
            rgauss(i,5)=Gaussout(3); %sig_xy
            rgauss(i,6)=Gaussout(4); %sig_z
    else
        rgauss(i,:) = NaN;
    end        
    
end

%%
fh = figure; 
height_offset = 0.1;
plot_height = (1 - 4*height_offset)/3;

surface_height = max(rgauss(:,3)*zvox/1000);

%I vs z plot
axes1 = axes('Parent',fh,'Position',[0.2,3*height_offset+2*plot_height,0.6,plot_height],'FontSize',12);
hold(axes1,'on');

plot(surface_height - rgauss(:,3)*zvox/1000,rgauss(:,4),...
    'MarkerSize',5,'Marker','o','LineStyle','none',...
    'Parent',axes1,'MarkerFaceColor',[0 0.45 0.74]);
rs = sortrows(rgauss,3);
plot(surface_height - rs(:,3)*zvox/1000,smooth(rs(:,3),rs(:,4),500,'rlowess'),...
    'LineWidth',2,...
    'Parent',axes1,'MarkerFaceColor',[1 0 0]);

ylim([0,1.5e6]);
xlim([-10,50]);
xlabel('Depth in microns from top of acquired volume');
ylabel('Spot Intensity');

%sigma_xy vs z plot
axes2 = axes('Parent',fh,'Position',[0.2,2*height_offset+plot_height,0.6,plot_height],'FontSize',12);
hold(axes2,'on');

plot(axes2,surface_height - rgauss(:,3)*zvox/1000,xvox*rgauss(:,5),...
    'MarkerSize',5,'Marker','o','LineStyle','none',...
    'MarkerFaceColor',[0.93 0.84 0.84],...
    'Color',[0.93 0.84 0.84]);

plot(surface_height - rs(:,3)*zvox/1000,xvox*smooth(rs(:,3),rs(:,5),500,'rlowess'),...
    'LineWidth',2,...
    'Parent',axes2,'MarkerFaceColor',[1 0 0]);
xlim([-10,50]);
xlabel('Depth in microns from top of acquired volume');
ylabel('Spot radial Size \sigma_{xy} (nm)');

%sigma_z vs z plot
axes3 = axes('Parent',fh,'Position',[0.2,height_offset,0.6,plot_height],'FontSize',12);
hold(axes3,'on');

plot(axes3,surface_height - rgauss(:,3)*zvox/1000,zvox*rgauss(:,6),...
    'MarkerSize',5,'Marker','o','LineStyle','none',...
    'MarkerFaceColor',[0.95 0.87 0.73],...
    'Color',[0.95 0.87 0.73]);

plot(surface_height - rs(:,3)*zvox/1000,zvox*smooth(rs(:,3),rs(:,6),500,'rlowess'),...
    'LineWidth',2,...
    'Parent',axes3,'MarkerFaceColor',[1 0 0]);
xlim([-10,50]);
xlabel('Depth in microns from top of acquired volume');
ylabel('Spot Axial Size \sigma_{z} (nm)');
