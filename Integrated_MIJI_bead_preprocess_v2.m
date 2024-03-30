%% Bead Pre process script
% -MIJI install
% -FIJI

% IMPORTANT all bead acquisition folders must have the word 'beads' in the folder name


% Also WD keeps changing away. Need to enforce it. or not uwe pwd anywhere

% add path to Fiji matlab scripts
addpath("C:\Users\Lionnet Lab\Documents\Fiji.app\scripts");
% download mij.jar and move to matlab jar folder
% make sure you have ij.jar 
% add java paths
javaaddpath("C:\Program Files\MATLAB\R2023b\java\jar\ij-1.54f.jar");
javaaddpath("C:\Program Files\MATLAB\R2023b\java\jar\mij.jar");
% run miji!
Miji();
MIJ.help %tells u what commands u can do
% which system used?


% add path to dependencies
% %% USER INPUT
% parentPath = 'D:\Finn\20240229_Staggered_PPARG_DNAFISH_Anlaysis';
% nChannels = 4; % number of channels in ur hyperstacks
% doFlatField =1;
% process_images = 1;
% 
% ref_channel = 2; % the channel you will map your data to. Use uManager channel value, eg 'C3-beads.loc3'
% 
% voxSize = [98, 98, 250]; % x,y,z of your voxel in nm
process_images = 1;
%% Search for fish img directories



all_FISH_dirs = {}; % this will contain images to split and process



subDirCtr = 0;
inImgDir = 0; % goes to 1 when we found the sub dir with images
% find folders in parentPath called FISHimgs

all_files = dir(fullfile(parentPath, '**/*.*'));

all_files = struct2table(all_files);

% get sub table with only images
img_idx = contains(all_files.name, 'tif');
img_table = all_files(img_idx, :);

% get unique dataset directory names
sampleDirs = unique(img_table.folder);

% keep only beads and remove if there is already a split_channels folder
bead_keep_idx = [];
for i = 1:numel(sampleDirs)

    if contains(string(sampleDirs{i}), 'bead', 'IgnoreCase',true) && ~contains(string(sampleDirs{i}), 'split_channels', 'IgnoreCase',true)
        disp('keeping only bead images');
        
     
        bead_keep_idx = [bead_keep_idx; i]; % get cond dir to report to user
    end


end

sampleDirs= sampleDirs(bead_keep_idx);

disp('');
disp(string(numel(sampleDirs)) + ' Unique image datasets were found:');
disp(string(sampleDirs));
disp('~~~~~~~~~');
disp('');

%%
% check if we already made FISH_img_Dirs and cfgs
skip_file_org = 0;
make_cfg = 1;
for i = 1:numel(sampleDirs)

    if contains(sampleDirs{i}, 'FISH_imgs')
        disp('FISH_imgs dirs are already present');
        skip_file_org = 1;
        
        all_FISH_dirs{i,1} = sampleDirs{i}; % save the FISH img dir
        sampleDirs{i} = fileparts(sampleDirs{i}); % get cond dir to report to user
    end


end

% see if we have made a cfg file yet for the beads
if sum(contains(all_files.name, '.ini') & contains(all_files.name, 'bead')) ~=0
    disp('cfg files already present');
    skip_cfg = 1;
end

% if we have already consolidated images in a subdir, make sure we define sample dirs
% by the parent diectory name
for i = 1:numel(sampleDirs)
    if contains(sampleDirs{i}, 'FISH_imgs')
        
        parent = fileparts(sampleDirs{i});
        sampleDirs{i} = parent;

    end

end

disp('');
disp(string(numel(sampleDirs)) + ' Unique image datasets were found:');
disp(string(sampleDirs));
disp('~~~~~~~~~');
disp('');

%%

% create a FISH img dir inside each unique sample dir
for i = 1:numel(sampleDirs)

    disp('');   
    disp('~~~~~~~~~');

    disp('Processing files in dir: ');
    disp(string(sampleDirs{i}));
    disp('');

    curFISHdir = fullfile(sampleDirs{i}, 'FISH_imgs');
    
    % make the fish imgs dir if it does not exist
    if ~exist(curFISHdir, 'dir')
        mkdir(curFISHdir);
    end
    
    % add FISH dir to list if not already present
    if ~skip_file_org
        all_FISH_dirs = vertcat(all_FISH_dirs, curFISHdir);
    end

    
    
    
    % get images for current sample
    if ~skip_file_org
        my_FOVpaths = img_table(strcmp(img_table.folder, sampleDirs{i}), 'name');
    else
        my_FOVpaths = img_table(strcmp(img_table.folder, curFISHdir), 'name');

    end

    
    
    % define root name for images
    [~,root_name] = fileparts(sampleDirs{i});
    root_name = string(root_name);
%     disp('Renaming files as: ' + string(root_name));
%     disp('Moving files to: ' + string(curFISHdir));
    disp('');   
    disp(string(numel(my_FOVpaths.name)) + ' FOVs in this dataset');

%     loop thru all images in current sample folder
    
% skip organizatin if already done
    if ~skip_file_org
        disp('Creating img dirs and renaming files');
        for j = 1:numel(my_FOVpaths.name)
            % give it a name with FOV index at the end
            new_name = root_name + '_' + string(j) + '.tif';
    
            old_path = fullfile(sampleDirs{i}, my_FOVpaths.name{j});
            new_path = fullfile(curFISHdir, new_name);
            
            % move all indexed FOVs to the FISH img folder
            movefile(old_path, new_path);
        end
    end

    if make_cfg == 1

        %   make a cfg file
        cfg_writer_for_pipe(sampleDirs{i}, root_name, nChannels, ref_channel, voxSize);
        disp('cfg file written to ' + string(fileparts(sampleDirs{i})));
    end

end

disp('~~~~~');
disp('Pre-processing Complete');
disp('~~~~~');
%% loop thru sample folders and split stacks


if process_images ==  1
    disp('');
    disp('~~~~~~~~~');
    disp('Splitting multi-channel images');
    disp('');
    
    
    % loop thru FISH img dirs one condition at a time
    for i = 1:numel(all_FISH_dirs)
    
        [curSampleDir,~] = fileparts(all_FISH_dirs{i});
    
        % make split channel dir
        split_ch_dir = fullfile(curSampleDir, 'split_channels');
        if ~exist (split_ch_dir, 'dir')
            mkdir(split_ch_dir)
        end
    
    % 
        % get all FOVS in cur fish dir FOVs
        cur_FOVs = dir(fullfile(all_FISH_dirs{i}, '*.tif'));
    
        cur_FOVs = struct2table(cur_FOVs);
        
        
        % loop thru and split
        for j = 1:numel(cur_FOVs(:,"name"))
    
            if numel(cur_FOVs(:,"name")) == 1
    
                my_FOV = cur_FOVs.name;
            
                my_FOVdir = cur_FOVs.folder;
        
                my_FOVpath = fullfile(string(my_FOVdir) , string(my_FOV));
            
            else
           
                my_FOV = cur_FOVs.name(j);
                
                my_FOVdir = cur_FOVs.folder(j);
        
                my_FOVpath = fullfile(string(my_FOVdir) , string(my_FOV));
            end
            %MIJI open my_FOVpath 
            
            disp('');
            
            disp('Processing')
            disp(my_FOVpath);
            
            % close all open windows
            MIJ.closeAllWindows();
            
            % open the image in FIJI
            MIJ.run('Open...', strcat('path=', my_FOVpath));
            
            % get the title
		    curTitle = char(MIJ.getCurrentTitle());
    
            MIJ.run('')
            
            % split channels
            MIJ.run('Split Channels');
            
            splitImg_list = MIJ.getListImages;
            
            nCh = numel(splitImg_list);
     
            for m=1:nCh
                
                
    
                
                MIJ.selectWindow(strcat('C', string(m), '-', curTitle));
                
                if doFlatField ==1 
                    MIJ.run('Pseudo flat field correction', 'blurring=50 hide stack');
                end
                
                saveTitle = char(MIJ.getCurrentTitle);
                MIJ.run('Save', strcat('Tiff..., path=[', fullfile(split_ch_dir, saveTitle), ']'));
                   
                
              
                
            end
        end
        
    
    end
end

MIJ.closeAllWindows();
MIJ.exit();

    




