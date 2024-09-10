%% FRIENDLY FISH Pre-processing
% Finn Clark | 6/24/24
% An integrated FISH pre procesing pipeline
% just press run and follow the GUI
% You must have a project folder set up as follows

% PROJECT DIR
% -FISH_condition1
% --img1.tif
% --img2.tif
% -FISH_condition2
% --img1.tif
% --img2.tif
% -my_project_beads (if using bead correction)
% --bead_img1.tif
% --bead_img2.tif

% Your bead folder must have the word 'bead' in it
% Your FISH conditions cannot have the word 'bead in it'
% all images are hypestacks saved by MM 
% dont worry about stack names, they will be renamed based on conditions


%% GUI Time
code_dir = fileparts(matlab.desktop.editor.getActiveFilename);
cd(code_dir)

% Generate a path string including all subfolders
pathWithSubfolders = genpath(code_dir);

% Add the path to MATLAB's search path
addpath(pathWithSubfolders);

% Display a message indicating that the path has been added
disp(['Added the following directory to the MATLAB path: ' code_dir]);

% See if we have correct FIJI path already saved in our txt file
FIJI_path_text = fullfile('FIJI_path.txt');
FIJI_path = fileread(FIJI_path_text);

FIJI_path = checkAndPromptForFIJIPath(FIJI_path, code_dir);

%~~~~RUN~~~~%
MIJ_installer


%% Assign Settings GUI
%~~~~RUN~~~~%
PreProcess_GUI
voxSize = [voxSizeXY, voxSizeXY, voxSizeZ];
% check for bad paths
checkFilePathForSpaces(project_dir);
checkFilePathForSpaces(code_dir);


%% settings

% mask_chan = 1; % specify uManager channel index to be used for masking (these will not be flatfield corrected)
do_maxZ = 1; % make max z projections of the mask_chans
% make_cell_masks = 1;
make_cfg = 1;

process_images = 1;
% doFlatField =1; 

% in pixels
% avg_cell_diameter= 150; 
% avg_cell_diameter= 70; b cells

cp_thresh = 0;


% You can use your own cell masks but they must be named according to
% the formula: 'conditionName_FOV.tif_ZMAX_MASK.tif
    % yes I realize this kinda sucks

%% RUN

DNA_FISH_Preprocess_Pipeline

