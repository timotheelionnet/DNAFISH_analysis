addpath('../fileSet');
addpath('../iniconfig');
addpath('../AIRLOCALIZE1_stable/AIRLOCALIZE1_4_export');
addpath('../AIRLOCALIZE1_stable/AIRLOCALIZE1_4_export/AIRLOCALIZE_1_4_subfunctions');
addpath('../salina');
addpath('../useful');
%addpath('../DiSPIM_PSF'); %This folder is not avialible on server and not
%used in this script
%% get the list of files
try
    configFileName;
catch
    % config file path
    configFileName = '20210111_DNA_FISH_crop_config.ini';
end

% create file list object
fs = fileSet;
fs.buildFileSetFromConfig(configFileName);

% reformat channel, replicate and FOV as numbers
convertNonNumeralStrings = 1; 
fs.toNumeric({'Channel','TechnicalReplicate','FOV'},convertNonNumeralStrings);

%% read parameters from config file relative to analysis into a parameters object
params = readParamsFromIniFile(configFileName);

%% loop through FISH Channels and perform AIRLOCALIZE
fs.addFileType('loc3',[]);
if params.channelDescription.useAirlocalizeGUI
    for i=1:numel(params.channelDescription.fishChannelList)
        curChannel = params.channelDescription.fishChannelList(i);
        curFileList = fs.getFileName({'Channel'},{curChannel},'fishImg','all');
        for j=1:numel(curFileList)
            disp('Processing file:')
            disp(curFileList{j})
             if j==1
                % run airlocalize with GUI the first time
                res = AIRLOCALIZE(curFileList{j});
                locpar = res.params;
                locpar.set_psf_manually = 0;
                locpar.set_thresh_manually = 0;
             else
                 % load current image
                 locpar.data = readTifStackWithImRead(curFileList{j});

                 % run airlocalize without GUI for the rest of the files
                 [res,~] = perform_detection_on_single_image_once_Dipankar2(locpar);
             end 

            loc3 = res.final_pix;
             
            % build the loc filename to enter in the fileSet list
            [d,f,~] = fileparts(curFileList{j});
            locFileName = fullfile(d,[f,'.loc3']);

            % find the conditions in the current file
            fl = fs.getFileName({},{},'fishImg','all');
            idx = find( ismember(fl,curFileList{j}) );
            curCond = {};
            for k = 1:numel(fs.conditions)
                curCond = [curCond,fs.fList.(fs.conditions{k})(idx)];
            end
            fs.setFileName(fs.conditions,curCond,'loc3',locFileName);
            
            save(locFileName,'loc3','-ascii');
        end
    end
else
    disp('GUI-less option not supported yet!');
end

%% loop through FISH Channels and plot dist hist
nDims = 3;

% compute a histogram of randomly localized spots for histogram normalization
% find image size
curChannel = params.channelDescription.fishChannelList(1);
curImgFileName = fs.getFileName({'Channel'},{curChannel},'fishImg','first');
curImgFileName = curImgFileName{1};
imSize = getImageSizeFromFile(curImgFileName); % get image size in pixels

% convert to nm
voxSize = [params.cropSpotsSettings.voxSize_dxy,...
            params.cropSpotsSettings.voxSize_dxy,...
            params.cropSpotsSettings.voxSize_dz];
imSize_nm = convert_loc_pix_to_nm(imSize,voxSize);

[nRand,randBinsCenters] = generate_pairDist_hist_randSpots(...
    params.cropSpotsSettings.nRandSpots,imSize_nm,...
    params.cropSpotsSettings.distMatrixHistogram_binSize);

% generate output folder for histogram figures
outFolderHist = fullfile(params.Output.outFolder,'pairDistHist');
if ~ exist(outFolderHist,'dir')
    mkdir(outFolderHist);
end

for i=1:numel(params.channelDescription.fishChannelList)
    curChannel = params.channelDescription.fishChannelList(i);
    curFileList = fs.getFileName({'Channel'},{curChannel},'loc3','all');
    for j=1:numel(curFileList) 
        
        disp('Processing file:')
        disp(curFileList{j})
            
        % load loc file
        loc = load(curFileList{j},'-ascii');
        
        % convert to nm (assumes 3D)
        loc3nm = convert_loc_pix_to_nm(loc,voxSize);
        
        % compute distance matrix
        [dr, ~] = compute_distMatrix_from_loc(loc3nm,[],nDims);
        
        % compute pair distance histogram
        [n,x] = hist(dr( ~isnan(dr) ), randBinsCenters);
        
        % normalize
        fh = figure('name','Normalized pair distances between spots');
        plot(x,n./nRand);
        xlabel('Distance in nm');
        ylabel('Number of spot pairs');
        
        % save histogram
        [~,f,~] = fileparts(curFileList{j});
        savefig(fh,fullfile(outFolderHist,[f,'_pairDistHist.fig']));
        close(fh);
    end
end

%% loop through FISH Channels and clean up double detections

fs.addFileType('cleanLoc3',[]);

outFolderCleanup = fullfile(params.Output.outFolder,'Cleaned');
if ~ exist(outFolderCleanup,'dir')
    mkdir(outFolderCleanup);
end

for i=1:numel(params.channelDescription.fishChannelList)
    curChannel = params.channelDescription.fishChannelList(i);
    curFileList = fs.getFileName({'Channel'},{curChannel},'loc3','all');
    for j=1:numel(curFileList) 
        
        disp('Processing file:')
        disp(curFileList{j})
        
        % load loc file
        loc = load(curFileList{j},'-ascii');
        
        % convert to nm (assumes 3D)
        loc3nm = convert_loc_pix_to_nm(loc,voxSize);
        
        % clean up spots
        cleanSpots = clean_up_double_detection(loc3nm, params.cropSpotsSettings.distThreshold, nDims);
        
        % output cleaned spot list
        [d,f,~] = fileparts(curFileList{j});
        locFileName = fullfile(d,[f,'_clean.loc3']);
        
        % output reconstructed image
        pos = convert_loc_pix_to_nm(cleanSpots,1./voxSize);
        save(locFileName,'pos','-ascii');

        imSize2 = [imSize(2),imSize(1),imSize(3)];
        spotSize = [params.cropSpotsSettings.spotSizeInImg,params.cropSpotsSettings.spotSizeInImg,params.cropSpotsSettings.spotSizeInImg];
        s = generate_gauss_spots_in_stack(pos,spotSize,imSize2);
        save_as_tiff(s,fullfile(outFolderCleanup,[f,'_clean.tiff']))
        
        % find the conditions in the current file
        fl = fs.getFileName({},{},'loc3','all');
        idx = find( ismember(fl,curFileList{j}) );
        curCond = {};
        for k = 1:numel(fs.conditions)
            curCond = [curCond,fs.fList.(fs.conditions{k})(idx)];
        end
        fs.setFileName(fs.conditions,curCond,'cleanLoc3',locFileName);
        
    end
end

%% Collect boxes and correct background 
[wx,wy,wz] = deal(params.cropSpotsSettings.boxSize_xy,params.cropSpotsSettings.boxSize_xy,params.cropSpotsSettings.boxSize_z);
thickness = params.cropSpotsSettings.boxBackgroundThickness;
fs.addFileType('crop',[]);

for i=1:numel(params.channelDescription.fishChannelList)
    curChannel = params.channelDescription.fishChannelList(i);
    curFileList = fs.getFileName({'Channel'},{curChannel},'fishImg','all');
    for j=1:numel(curFileList) 
        disp('Processing file:')
        disp(curFileList{j})
        
        %load img and loc3
        ims = readTifStackWithImRead(curFileList{j});
        [d,f,~] = fileparts(curFileList{j});
        loc = load(fullfile(d, [f, '_clean.loc3']),'-ascii');
        
        %create output folder        
        outFolderCrop = fullfile(params.Output.outFolder,[f,'_crop']);
        if ~ exist(outFolderCrop,'dir')
            mkdir(outFolderCrop);
        end
        fl = fs.getFileName({},{},'fishImg','all');
        idx = find( ismember(fl,curFileList{j}) );
        curCond = {};
        for k = 1:numel(fs.conditions)
            curCond = [curCond,fs.fList.(fs.conditions{k})(idx)];
        end
        fs.setFileName(fs.conditions,curCond,'crop',outFolderCrop);
        
        %crop and output boxes
        boxSize = [params.cropSpotsSettings.boxSize_xy, params.cropSpotsSettings.boxSize_xy, params.cropSpotsSettings.boxSize_z];
        if params.cropSpotsSettings.testing
            crop_box(ims, loc, outFolderCrop, boxSize, params.cropSpotsSettings.boxBackgroundThickness);
        else
            crop_box_hotpixel(ims, loc, outFolderCrop, boxSize, params.cropSpotsSettings.boxBackgroundThickness);
        end
    end
end
