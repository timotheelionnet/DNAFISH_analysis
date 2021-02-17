%% Collect all .loc files from one analysis,
%loads them, removes negative spots, generate a reconstructed z-stack with
%the resulting spots (for visual inspection) and compiles the number of spots in one table.

%plots intensity distributions for each condition

%plots the copy number per cell as a plotSpread for each condition

%IMPORTANT:
%The files analyzed here correspond to segmented cells in 5x Z-compressed
%datasets.
%1) datasets were compressed in Z by a factor of 5 by (?) performing a max
    %projection of bloacks of 5 slices.
%2) cells were manually segmented on the compressed dataset using TrackEM,
%3) each segmented cell was saved as a separate (Z-compressed) z-stack,
    %with the regions outside of the segmented ROI set to zero
%4) individual segmented stacks (one per cell) were analyzed using AIRLOCALIZE with the following
%parameters:
    % sigma_xy:    0.65304 
    % sigma_z:    0.21341 
    % dx:    64 
    % dy:    64 
    % dz:    200 
    % cutsize:    3 
    % tol:    0.01 
    % thickness:    1 
    % thresh.level:    3000 
    % thresh.units:    absolute 
    % thresh.sd:    139048.8964 
    % max_spots:    200000 
    % ROIsize:    2 
    % method:    could not read value 
    % fit:    could not read value 
    % maxcount:    100 
    % data_fname:     
    % data_stackname:    could not read value 
    % data_dirname:    could not read value 
    % save_dirname:    could not read value 
    % numdim:    3 
    % filter.nlo:    2.15 
    % filter.nhi:    1 
    % type:    could not read value 
    % set_psf_manually:    1 
    % set_thresh_manually:    1 
    % bg_mode:    could not read value 
    % mode:    could not read value 
    % AIRLOCALIZE_version:    1 
    % output_spot_image:    0 
    % dense_spots:    0 
    % dirmode:    could not read value 
    % select_ROI:    0 
    % nx:    66 
    % ny:    199 
    % nz:    24 
    % curfname:    could not read value 

%% set matlab path to include necessary subfunctions
%with platform dependent path names

%Launch_from = 'cluster';
Launch_from = 'desktop';

switch Launch_from
    case 'cluster'
        data_rootname = '/groups/singcellbio/';
        loc_rootname = '/groups/transimcon/home/';
    case 'desktop'
        data_rootname = '/Volumes/';
        loc_rootname = '/Volumes/';
end

addpath([loc_rootname,'lionnett/MatlabCode/plotSpread']);

%intensity threshold to discard spots
signal_threshold = 0;

%% Large ZT2
%collect file names
folder_large_tif_zt2 = '/Volumes/transimcon/Salina/Data/Timeless SM/timeless ZT2 BB+SIM/Timeless ZT2 Large LNvs/';

tif_flist_ZT2_large = get_clean_file_list(folder_large_tif_zt2,{'C1','tif'}, {'MAX','recon.tif'},0,0);

for i=1:numel(tif_flist_ZT2_large)
    loc_flist_ZT2_large{i} = strrep(tif_flist_ZT2_large{i},'tif','loc3');
end

% loop through ZT 2 large data to load detections, generate and save reconstructed z-stack
spots_ZT2_large = [];
spots_sig_ZT2_large = [];
n_spots_ZT2_large = [];
n_spots_with_bg_ZT2_large = [];
for i=1:numel(loc_flist_ZT2_large)
    
    disp([num2str(i),'/',num2str(numel(loc_flist_ZT2_large)),'; ',loc_flist_ZT2_large{i}]);
    spots_ZT2_large{i} = dlmread(loc_flist_ZT2_large{i});
    stack = timtiffread(tif_flist_ZT2_large{i});
    stacksize = size(stack);
    clear stack;
    
    spots_sig_ZT2_large{i} = spots_ZT2_large{i}(spots_ZT2_large{i}(:,4)>signal_threshold,:);
    sp = spots_sig_ZT2_large{i};
    sp(:,4) = sp(:,4)/10;
    stack_sig = generate_gauss_spots_in_stack(sp,[1,1,1],stacksize);

    [f,n,e] = fileparts(tif_flist_ZT2_large{i});

    save_as_tiff(stack_sig, fullfile(f,[n,'_recon.tif']));
    n_spots_ZT2_large(i) = size(spots_sig_ZT2_large{i},1);
    n_spots_with_bg_ZT2_large(i) = size(spots_ZT2_large{i},1);
end

%% combine all spots from the ZT2 large condition into one table (rather than a cell array)
spots_sig_ZT2_large_combined= [];
for i=1:numel(spots_sig_ZT2_large)
    spots_sig_ZT2_large_combined = [spots_sig_ZT2_large_combined; [spots_sig_ZT2_large{i},i*ones(size(spots_sig_ZT2_large{i},1),1)]];
    
end

[nInt_ZT2_large,xInt_ZT2_large] = hist(spots_sig_ZT2_large_combined(:,4),0:1000:max(spots_sig_ZT2_large_combined(:,4)));
f_int_hist = figure('Name','Spot Intensity Distribution Timeless (5x z-Compressed data)');
hold;
a_Int_hist = gca;
plot(a_Int_hist,xInt_ZT2_large,nInt_ZT2_large/sum(nInt_ZT2_large),...
    'DisplayName','ZT2 large, Spot Intensity Distribution');
plot(a_Int_hist,xInt_ZT2_large,cumsum(nInt_ZT2_large/sum(nInt_ZT2_large)),...
    'DisplayName','ZT2 large, Spot Intensity CDF');
xlabel('Spot Intensity (counts)');
ylabel('Probability / Cumulated Probability');

%% Small ZT2
%collect file names
folder_tif_ZT2_small = '/Volumes/transimcon/Salina/Data/Timeless SM/timeless ZT2 BB+SIM/Timeless ZT2 Small LNvs/';

tif_flist_ZT2_small = get_clean_file_list(folder_tif_ZT2_small,{'C1','tif'}, {'MAX','recon.tif'},0,0);

for i=1:numel(tif_flist_ZT2_small)
    loc_flist_ZT2_small{i} = strrep(tif_flist_ZT2_small{i},'tif','loc3');
end

% loop through ZT 2 small data to load detections, generate and save reconstructed z-stack
spots_ZT2_small = [];
spots_sig_ZT2_small = [];
n_spots_ZT2_small = [];
n_spots_with_bg_ZT2_small = [];
for i=1:numel(loc_flist_ZT2_small)
    
    disp([num2str(i),'/',num2str(numel(loc_flist_ZT2_small)),'; ',loc_flist_ZT2_small{i}]);
    spots_ZT2_small{i} = dlmread(loc_flist_ZT2_small{i});
    stack = timtiffread(tif_flist_ZT2_small{i});
    stacksize = size(stack);
    clear stack;
    
    spots_sig_ZT2_small{i} = spots_ZT2_small{i}(spots_ZT2_small{i}(:,4)>signal_threshold,:);
    sp = spots_sig_ZT2_small{i};
    sp(:,4) = sp(:,4)/10;
    stack_sig = generate_gauss_spots_in_stack(sp,[1,1,1],stacksize);

    [f,n,e] = fileparts(tif_flist_ZT2_small{i});

    save_as_tiff(stack_sig, fullfile(f,[n,'_recon.tif']));
    n_spots_ZT2_small(i) = size(spots_sig_ZT2_small{i},1);
    n_spots_with_bg_ZT2_small(i) = size(spots_ZT2_small{i},1);
end

%% combine all spots from the ZT2 small condition into one table (rather than a cell array)
spots_sig_ZT2_small_combined= [];
for i=1:numel(spots_sig_ZT2_small)
    spots_sig_ZT2_small_combined = [spots_sig_ZT2_small_combined; [spots_sig_ZT2_small{i},i*ones(size(spots_sig_ZT2_small{i},1),1)]];
    
end
[nInt_zt2_small,xInt_zt2_small] = hist(spots_sig_ZT2_small_combined(:,4),0:1000:max(spots_sig_ZT2_small_combined(:,4)));

plot(a_Int_hist,xInt_zt2_small,nInt_zt2_small/sum(nInt_zt2_small),...
    'DisplayName','ZT2 small, Spot Intensity Distribution');
plot(a_Int_hist,xInt_zt2_small,cumsum(nInt_zt2_small/sum(nInt_zt2_small)),...
    'DisplayName','ZT2 small, Spot Intensity CDF');


%% ZT14 small
%collect file names
folder_tif_ZT14_small = '/Volumes/transimcon/Salina/Data/Timeless SM/timeless ZT2 BB+SIM/Timeless ZT14 Small LNvs/';

tif_flist_ZT14_small = get_clean_file_list(folder_tif_ZT14_small,{'C1','tif'}, {'MAX','recon.tif'},0,0);

for i=1:numel(tif_flist_ZT14_small)
    loc_flist_ZT14_small{i} = strrep(tif_flist_ZT14_small{i},'tif','loc3');
end

% loop through ZT 14 small data to load detections, generate and save reconstructed z-stack
spots_ZT14_small = [];
spots_sig_ZT14_small = [];
n_spots_ZT14_small = [];
n_spots_with_bg_ZT14_small = [];
for i=1:numel(loc_flist_ZT14_small)
    
    disp([num2str(i),'/',num2str(numel(loc_flist_ZT14_small)),'; ',loc_flist_ZT14_small{i}]);
    spots_ZT14_small{i} = dlmread(loc_flist_ZT14_small{i});
    stack = timtiffread(tif_flist_ZT14_small{i});
    stacksize = size(stack);
    clear stack;
    
    spots_sig_ZT14_small{i} = spots_ZT14_small{i}(spots_ZT14_small{i}(:,4)>signal_threshold,:);
    sp = spots_sig_ZT14_small{i};
    sp(:,4) = sp(:,4)/10;
    stack_sig = generate_gauss_spots_in_stack(sp,[1,1,1],stacksize);

    [f,n,e] = fileparts(tif_flist_ZT14_small{i});

    save_as_tiff(stack_sig, fullfile(f,[n,'_recon.tif']));
    n_spots_ZT14_small(i) = size(spots_sig_ZT14_small{i},1);
    n_spots_with_bg_ZT14_small(i) = size(spots_ZT14_small{i},1);
end

%% combine all spots from the ZT14 small condition into one table (rather than a cell array)
spots_sig_ZT14_small_combined= [];
for i=1:numel(spots_sig_ZT14_small)
    spots_sig_ZT14_small_combined = [spots_sig_ZT14_small_combined; [spots_sig_ZT14_small{i},i*ones(size(spots_sig_ZT14_small{i},1),1)]];
    
end
[nInt_zt14_small,xInt_zt14_small] = hist(spots_sig_ZT14_small_combined(:,4),0:1000:max(spots_sig_ZT14_small_combined(:,4)));

plot(a_Int_hist,xInt_zt14_small,nInt_zt14_small/sum(nInt_zt14_small),...
    'DisplayName','ZT14 small, Spot Intensity Distribution');
plot(a_Int_hist,xInt_zt14_small,cumsum(nInt_zt14_small/sum(nInt_zt14_small)),...
    'DisplayName','ZT14 small, Spot Intensity CDF');

%% ZT14 large
%collect file names
folder_tif_ZT14_large = '/Volumes/transimcon/Salina/Data/Timeless SM/timeless ZT2 BB+SIM/Timeless ZT14 Large LNvs/';

tif_flist_ZT14_large = get_clean_file_list(folder_tif_ZT14_large,{'C1','tif'}, {'MAX','recon.tif'},0,0);

for i=1:numel(tif_flist_ZT14_large)
    loc_flist_ZT14_large{i} = strrep(tif_flist_ZT14_large{i},'tif','loc3');
end

% loop through ZT 14 large data to load detections, generate and save reconstructed z-stack
spots_ZT14_large = [];
spots_sig_ZT14_large = [];
n_spots_ZT14_large = [];
n_spots_with_bg_ZT14_large = [];
for i=1:numel(loc_flist_ZT14_large)
    
    disp([num2str(i),'/',num2str(numel(loc_flist_ZT14_large)),'; ',loc_flist_ZT14_large{i}]);
    spots_ZT14_large{i} = dlmread(loc_flist_ZT14_large{i});
    stack = timtiffread(tif_flist_ZT14_large{i});
    stacksize = size(stack);
    clear stack;
    
    spots_sig_ZT14_large{i} = spots_ZT14_large{i}(spots_ZT14_large{i}(:,4)>signal_threshold,:);
    sp = spots_sig_ZT14_large{i};
    sp(:,4) = sp(:,4)/10;
    stack_sig = generate_gauss_spots_in_stack(sp,[1,1,1],stacksize);

    [f,n,e] = fileparts(tif_flist_ZT14_large{i});

    save_as_tiff(stack_sig, fullfile(f,[n,'_recon.tif']));
    n_spots_ZT14_large(i) = size(spots_sig_ZT14_large{i},1);
    n_spots_with_bg_ZT14_large(i) = size(spots_ZT14_large{i},1);
end

%% combine all spots from the ZT14 large condition into one table (rather than a cell array)
spots_sig_ZT14_large_combined= [];
for i=1:numel(spots_sig_ZT14_large)
    spots_sig_ZT14_large_combined = [spots_sig_ZT14_large_combined; [spots_sig_ZT14_large{i},i*ones(size(spots_sig_ZT14_large{i},1),1)]];
    
end
[nInt_zt14_large,xInt_zt14_large] = hist(spots_sig_ZT14_large_combined(:,4),0:1000:max(spots_sig_ZT14_large_combined(:,4)));
plot(a_Int_hist,xInt_zt14_large,nInt_zt14_large/sum(nInt_zt14_large),...
    'DisplayName','ZT14 large, Spot Intensity Distribution');
plot(a_Int_hist,xInt_zt14_large,cumsum(nInt_zt14_large/sum(nInt_zt14_large)),...
    'DisplayName','ZT14 large, Spot Intensity CDF');

%% plotspread - plot Timeless copy number per cell as a plotspread with one cloud per condition
%this includes all data points detected (except those with negative intensity)

fp = figure('Name','Timeless Copy Numbers (Inclusive Spot Intensity Threshold)');

plotSpread({n_spots_ZT2_small;n_spots_ZT2_large;n_spots_ZT14_small;n_spots_ZT14_large},...
    'distributionColors',{[1,0,0];[0,0,1];[1,0,0];[0,0,1];},...
    'distributionMarkers',{'o','o','o','o'},...
    'xValues',[1,2,3,4],...
    'spreadWidth',0.5,...
    'xNames',{'ZT2 sLNvs','ZT2 lLNvs','ZT14 sLNvs','ZT14 lLNvs'});

errorbar([1.35,2.25,3.25,4.25],...
[mean(n_spots_ZT2_small),mean(n_spots_ZT2_large),mean(n_spots_ZT14_small),mean(n_spots_ZT14_large)],...
[std(n_spots_ZT2_small),std(n_spots_ZT2_large),std(n_spots_ZT14_small),std(n_spots_ZT14_large)],'o');
ylabel('Timeless copy number per cell');

%% 'Conservative' plotspread - plot Timeless copy number per cell as a plotspread with one cloud per condition
%this includes data points detected with intensities greater than 50,000
    % sanity check because visual inspection reveals that there might be some false negatives on the
    % spots below 50000 counts. Turns out the results are not substantially
    % changed.
signal_threshold_conservative = 50000;

n_spots_ZT2_small_conservative = [];
for i=1:numel(spots_sig_ZT2_small)
    n_spots_ZT2_small_conservative(i) = sum(spots_sig_ZT2_small{i}(:,4)>signal_threshold_conservative);
end
n_spots_ZT2_large_conservative = [];
for i=1:numel(spots_sig_ZT2_large)
    n_spots_ZT2_large_conservative(i) = sum(spots_sig_ZT2_large{i}(:,4)>signal_threshold_conservative);
end
n_spots_ZT14_small_conservative = [];
for i=1:numel(spots_sig_ZT14_small)
    n_spots_ZT14_small_conservative(i) = sum(spots_sig_ZT14_small{i}(:,4)>signal_threshold_conservative);
end
n_spots_ZT14_large_conservative = [];
for i=1:numel(spots_sig_ZT14_large)
    n_spots_ZT14_large_conservative(i) = sum(spots_sig_ZT14_large{i}(:,4)>signal_threshold_conservative);
end

fp = figure('Name','Timeless Copy Numbers (Conservative Spot Intensity Threshold)');
plotSpread({n_spots_ZT2_small_conservative;...
    n_spots_ZT2_large_conservative;...
    n_spots_ZT14_small_conservative;...
    n_spots_ZT14_large_conservative...
    },...
    'distributionColors',{[1,0,0];[0,0,1];[1,0,0];[0,0,1];},...
    'distributionMarkers',{'o','o','o','o'},...
    'xValues',[1,2,3,4],...
    'spreadWidth',0.5,...
    'xNames',{'ZT2 sLNvs','ZT2 lLNvs','ZT14 sLNvs','ZT14 lLNvs'});

errorbar([1.35,2.25,3.25,4.25],...
[mean(n_spots_ZT2_small_conservative),mean(n_spots_ZT2_large_conservative),mean(n_spots_ZT14_small_conservative),mean(n_spots_ZT14_large_conservative)],...
[std(n_spots_ZT2_small_conservative),std(n_spots_ZT2_large_conservative),std(n_spots_ZT14_small_conservative),std(n_spots_ZT14_large_conservative)],'o');
ylabel('Timeless copy number per cell');

%plot mean vs STD
figure('Name','Mean vs STD Timeless');
plot([mean(n_spots_ZT2_small_conservative),mean(n_spots_ZT2_large_conservative),mean(n_spots_ZT14_small_conservative),mean(n_spots_ZT14_large_conservative)],...
    [std(n_spots_ZT2_small_conservative),std(n_spots_ZT2_large_conservative),std(n_spots_ZT14_small_conservative),std(n_spots_ZT14_large_conservative)].^2);
xlabel('Mean number of Timeless mRNA /cell');
ylabel('STD number of Timeless mRNA /cell');

%clear plot and axes handles
clear fp a_Int_hist f_int_hist

