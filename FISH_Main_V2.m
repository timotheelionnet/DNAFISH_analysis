%% DNA FISH ANALYSIS PIPELINE
% Lionnet Lab, 2024
% Input = DNA FISH hyperstacks split by channel
% Output = Spatial analysis

% Instructions: Fill out USER INPUT and then run one chunk at a time,
% adjusting settings as needed

% Currently can analyze one condition at a time

%% USER INPUT
% paste path to the fodler containing your beads (optional) and one folder
% for each FISH condition 
project_dir= 'H:\FC_DNA_FISH_PIPELINE_INSTALL_20240624\20240624_User_Friendly_FISH_Pipeline\demo';

%paste path to the ini file for the condition you would like to analyze
config_fname = "H:\FC_DNA_FISH_PIPELINE_INSTALL_20240624\20240624_User_Friendly_FISH_Pipeline\demo\MSC_pparg.ini";

% Load params
load_params_from_cfg_for_pipe

% Whether to use bead correction of chromatic aberration
perform_bead_correction = 1;
    % 1 (if you have already run bead pre-processing and analysis
    % 0 (you just want to analyze your raw FISH images


%% apply channel key
% Required files
% - channelKey.txt file--script first searches your sample dir, then the parent
% dir
% - fits.mat produced from BEAD analysis--script first searches your sample dir, then the parent
% dir

[key_found, key_path] = check_for_channel_key(project_dir);

% load key
key_tab = readtable(key_path);

% check for errors in key
verify_channel_key(key_tab)

% interpret channel key and write a cfg file
[~,parent] = fileparts(project_dir);
refChannel = key_tab.uM_FISH_ch( key_tab.isReference == 1);
bead_ref_ch = key_tab.uM_Bead_ch( key_tab.isReference == 1);
bead_ch_vals = key_tab.uM_Bead_ch( key_tab.uM_Bead_ch ~=0);


key_tab 

%% Organize bead correction functions

if perform_bead_correction == 1

    bead_dir = find_bead_dir(project_dir);
    %~~~RUN~~~%
    organize_bead_buddy_functions_for_pipe

    %~~~Run~~~%
    % apply_channel_key_V2
end




%% localize FISH images

%~~~User Settings~~~%

% 1) Choose whether to run AIRLOCALIZE (must run it first time analyzing)
runAIRLOCALIZE = 1; % If locpar are already loaded in your workspace: True = perform localization using existing locpar. False = skip. 

% 2) Option to save a tiff of reconstructed spots (very slow)
gen_gauss_spots = false;    
clean_double_detections = 0; % looks to see if spots in one channel are closer than our matching threshold. Keeps only the brighter spot

%~~~Run~~~%
localize_FISH_spots_for_pipe_V2;

%% apply reigstration functions to your data

if perform_bead_correction == 1

    % if you want to crop your data to a specifc ROI 
    crop = 0;
    roi_x = [200,1000];
    roi_y = [200,1000];
    roi = [roi_x; roi_y];
    
    
    %~~~Run~~~%
    register_Loc3_for_pipe_V2

end


%% sort spots into masks 
% required inputs
% - mask files and paths loaded into fs
% - loc3 files loaded by fileset 
% -- this could be direct from a cfg path that points to the reg loc3
% -- or this could be using the fs created by the previous step
rmv_bkrnd_spots = 1;

if perform_bead_correction == 1
    file2sort = 'loc_REG'; 
else
    file2sort = 'loc';
end

% set a number of spots that must be in the cell (need to validate this)
spotThresh = [1,4]; % at least min at most  max spots per channel
apply_spot_thresh = 1; 

%~~~Run~~~%
sort_loci_into_masks_for_pipe_V2


%% Spatial analysis and inter-channel spot matching

if perform_bead_correction == 1
    file2match = 'loc_REG_sorted'; % again change this to 'loc_sorted' if you did not do registration
else
    file2match = 'loc_sorted';
end

subtractAvgOffset = 0;  
    % set to 0 for raw coordinates (use this if you applied bead correction); 
    % set to 1 for channel-corrected coordinates
    % set to 2 for pair-corrected distances so that all median offsets
    % are strictly zero.

matchDistThreshold = 1000; %(nm)
    % spots must be within this distance to be called a match
    
mutualNN_only = 1;
    % --set to 0 to include all spot pairs within the threshold (allows spots
    % to have multiple match partners
    % --set to 1 to keep matches between mutual nearest neighbors (more
    % stringent and generally the way to go)


do_cdf_crxn = 0;
    % models the false posivite regime of your CDF and subtracts it from
    % your data. Corrected CDF table will be saved separately




trans_match_thresh = 5000; %(nm_)
    % matches beyond this distance are deemed false positves (physically
    % impossible)

%~~~Filter by connectivity~~~%
% these filters will be saved in a seprate results table in
% condition/res/matched_spots called connectivity_gated...etc
apply_spot_thresh = 1;
minSpotsPerGroup = 3;
maxSpotsPerGroup = 7;
minChannelsPerGroup = 3;
maxChannelsPerGroup = 4;
strictlyOneSpotPerChannel = 0;


%~~~Run~~~%

match_FISH_spots_V3


