%% cfg writer
% Did you move youre files around? Already processed images and performed
% detection?
% You can run this script to update your cfg files to match new file paths
parentPath = 'H:\Finn_DriveH\20240226_PPARG_Staggered_Timecourse\all_conds';
voxSize = [98, 98, 250]; % x,y,z of your voxel in nm
FISH_channels = [2,3,4,5];

make_cfg = 1;

sampleDirs = {}; % this will contain images to split and process

filePath = matlab.desktop.editor.getActiveFilename;
[code_dir,f,~] = fileparts(filePath);

cd(code_dir);

% Get a list of all items (files and directories) in the parent directory
allItems = dir(parentPath);

% Loop through all items and check if they are directories
for i = 1:length(allItems)
    if allItems(i).isdir && ~strcmp(allItems(i).name, '.') && ~strcmp(allItems(i).name, '..')
        % Exclude '.' and '..' directories
        % Check if the directory is one level deep
        if length(strfind(allItems(i).folder, filesep)) == length(strfind(parentPath, filesep))
            sampleDirs = [sampleDirs; fullfile(parentPath, allItems(i).name)];
        end
    end
end

% Display the list of directories
disp('List of immediate subdirectories:');
disp(sampleDirs);

bead_rmv_idx = [];
for i = 1:numel(sampleDirs)

    if contains(string(sampleDirs{i}), 'bead', 'IgnoreCase',true) || contains(string(sampleDirs{i}), 'split_channels', 'IgnoreCase',true)
        disp('removing bead dir');
        
     
        bead_rmv_idx = [bead_rmv_idx; i]; % get cond dir to report to user
    end


end

sampleDirs(bead_rmv_idx) = [];

disp('');
disp(string(numel(sampleDirs)) + ' Unique image datasets were found:');
disp(string(sampleDirs));
disp('~~~~~~~~~');
disp('');

for i = 1:numel(sampleDirs)
    % write a cfg file
    if make_cfg == 1

        %   make a cfg file
        [~,root_name] = fileparts(sampleDirs{i});
        FISH_cfg_writer_for_pipe(sampleDirs{i}, root_name, FISH_channels, voxSize)
        disp('cfg file written to ' + string(fileparts(sampleDirs{i})));
    end
end
