%Load_Trackmate_spots_and_save_stack

%% your data files here
data_xls = '/Users/lionnett/Documents/Salina Long/2016_08_24/segres/All Spots statistics 1px.xlsx';
data_xls = '/Users/lionnett/Documents/Salina Long/Mi149_data/Mi1_256/C4-Mi1 Sample6 R High Res_Linear unmixing-2_trackmate1.xlsx';
%should be xlsx format
%first sheet only
%data columns used:
%'POSITION_X'
%'POSITION_Y'
%'POSITION_Z'

%size of the stack in pixels
%use same axis order as in ImageJ
imsize = [850,850,159];
imsize = [256,256,159];

%this is the folder where the filtered data is saved
save_path = '/Users/lionnett/Documents/Salina Long/2016_08_24/segres/testBPfilt/';

%% dataset info


%size of the voxel in FIJI units
xpixsize = 1;
ypixsize = 1;
zpixsize = 1;

%size in pixels of cell representation in output visualizations [x,y,z]
cell_blob_size = [1,1,1];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

addpath([loc_rootname,'lionnett/MatlabCode/salina']);

%% load xls file
[num,txt,~] = xlsread(data_xls);

%% extract positions and intensities
header_strings = txt(1,:);

%locate column holding X,Y,Z positions, Intensities
X_Idx = find(ismember(header_strings,'POSITION_X'));
Y_Idx = find(ismember(header_strings,'POSITION_Y'));
Z_Idx = find(ismember(header_strings,'POSITION_Z'));

cell_pos_FIJIunits = num(:,[X_Idx,Y_Idx,Z_Idx]);
cell_pos_pixels = [cell_pos_FIJIunits(:,1)/xpixsize + 1, cell_pos_FIJIunits(:,2)/ypixsize + 1, cell_pos_FIJIunits(:,3)/zpixsize+1];

%% save data
stack_cells = generate_gauss_spots_in_stack(cell_pos_pixels,cell_blob_size,imsize);

%flip x/y to get the same conventions as in ImageJ
stack_cells = permute(stack_cells,[2,1,3]);

save_as_tiff(stack_cells, fullfile(save_path,'all_segmented_cells.tif' ));
