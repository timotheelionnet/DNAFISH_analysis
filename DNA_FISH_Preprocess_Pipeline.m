%% DNA FISH PRE-PROCESSING PIPELINE (run one chunk at a time)
% requires
% -MIJI install
% -FIJI with biovoxxel plug in
% -Add Ons: Computer Vision, Deep Learning, GPU Coder, Medical Imaging,
% Medical Image Interface for Cellpose

% ~~~~~~~Before starting, make sure that:~~~~~~~
% --You have added the pipeline directory to your path with all subfolders
% --all FOVs per condition are in the same folder 
% ----'beads' is not in FISH data folder name

% --your bead stacks are in a folder with the word 'beads' in it


%% USER INPUT: Project directory
% give the path to the project folder containing bead and/or FISH image folders
parentPath = 'G:\Finn\20240327_FLUOR_SWAP_PPARG\swapped_fluor_set';

%% USER INPUT: Pre-process Bead images

nChannels = 4; % number of channels in ur bead hyperstacks
doFlatField =1; % apply flat field correction to your iamges

ref_channel = 2; % the channel you will map your data to. Use uManager channel value in your beads, eg 'C3-beads.loc3'

voxSize = [98, 98, 250]; % x,y,z of your voxel in nm
%%
%~~~~RUN~~~~~~%
Integrated_MIJI_bead_preprocess_v2

%%  USER INPUT: Pre-process FISH images

mask_chan = 1; % specify uManager channel index in FISH images to be used for segmentation
do_maxZ = 1; % save max z projections of the mask_chans?
make_cell_masks = 1; % make cell masks with cellpose? Requires do_maxZ =1
make_cfg = 1; % save a config file for this sample? (Required for analysis)

process_images = 1; %do you want to actually process the images? Set 0 to only organize files
doFlatField =1; % apply flat field correction (I always use this)

% For cell segmentation
avg_cell_diameter= 150; % avg nuclear diameter in pixels
    % 70 pix for B-cells at 100x, 150 pix for MSCs at 60x Nikon

    % You can use your own cell masks but they must be named according to
    % the formula: 'conditionName_FOV.tif_ZMAX_MASK.tif
        % yes I realize this kinda sucks
%%
%~~~~RUN~~~~~~%
Integrated_MIJI_preProcess_v2


%% USER INPUT Bead analysis (must wait until pre-processing is done

% Above Preprocessing will generate a .ini file in your project directory for
% beads. Paste that path here and run chunk
beadConfigPath = "G:\Finn\20240327_FLUOR_SWAP_PPARG\beads.ini";

% output will be saved in beadfolder/res as a beadsName_fits.mat.  

%~~~RUN~~~%
Bead_Analysis_Pipeline_V6
