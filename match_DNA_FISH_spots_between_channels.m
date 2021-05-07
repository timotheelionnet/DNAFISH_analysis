addpath('../fileSet');
addpath('../iniconfig');

%% get list of files
try
    configFileName;
catch
    % config file path
    configFileName = '20210201_DNA_FISH_match_loci_config.ini';
end

% create file list object
fs = fileSet;
fs.buildFileSetFromConfig(configFileName);

% reformat channel, replicate and FOV as numbers
convertNonNumeralStrings = 1; 
fs.toNumeric({'Channel','TechnicalReplicate','FOV'},convertNonNumeralStrings);

%read parameters from config file relative to analysis into a parameters object
params = readParamsFromIniFile(configFileName);

%% loop through fish channels and match loci
if ~ exist(fullfile(params.Output.outFolder,'loci_match'),'dir')
    mkdir(fullfile(params.Output.outFolder,'loci_match'));
end
        
compute_hist = 0;
nDims = 3;
voxSize = [params.MatchLociSettings.voxSize_dxy,...
    params.MatchLociSettings.voxSize_dxy,...
    params.MatchLociSettings.voxSize_dz];

if compute_hist
    % compute a histogram of randomly localized spots for histogram normalization
    % find image size
    curChannel = params.channelDescription.fishChannelList(1);
    curImgFileName = fs.getFileName({'Channel'},{curChannel},'fishImg','first');
    curImgFileName = curImgFileName{1};
    imSize = getImageSizeFromFile(curImgFileName); % get image size in pixels

    % convert to nm
     imSize_nm = convert_loc_pix_to_nm(imSize,voxSize);

    [nRand,randBinsCenters] = generate_pairDist_hist_randSpots(...
        params.MatchLociSettings.nRandSpots,imSize_nm,...
        params.MatchLociSettings.distMatrixHistogram_binSize);

    % generate output folder for histogram figures
    outFolderHist = fullfile(params.Output.outFolder,'pairDistHist');
    if ~ exist(outFolderHist,'dir')
        mkdir(outFolderHist);
    end
end

fs.addFileType('loci_match',[]);

C2 = setdiff(params.channelDescription.fishChannelList, params.channelDescription.ReferChannel);

% loop through ReferChannel and assign cell ID to each loci
curChannel = params.channelDescription.ReferChannel;
curFileList = fs.getFileName({'Channel'},{curChannel},'fishImg','all');
DAPI = params.channelDescription.DAPI;

for j=1:numel(curFileList)
    % load loc file
    [~,f,~] = fileparts(curFileList{j});
    disp(['Assigning cell ID for ', f]);
    cropFolder1 = fullfile(params.InputFolders.CropFolder,[f,'_crop']);
    loc1 = load(fullfile(cropFolder1,'location.loc3'),'-ascii');
    if isempty(loc1)
        disp('Loc3 file is empty!')
        loc1 = [];
        save(fullfile(cropFolder1,'location_mask_sorted.loc3'),'loc1','-ascii');
        continue
    end
    %load cell mask
    [~,f,~] = fileparts(curFileList{j});
    MaskFile = fullfile(params.InputFolders.CellMaskFolder,strcat('C', string(DAPI),f(3:length(f)),'max_Mask.tiff'));
    Mask = readTifStackWithImRead(MaskFile);
    %sort loci into cell masks   
    loc1 = sort_loci_into_cell_masks(loc1,Mask);
    save(fullfile(cropFolder1,'location_mask_sorted.loc3'),'loc1','-ascii');
end

% loop through fish channels and match loci
for i=1:length(C2)
    curChannel = C2(i);
    curFileList = fs.getFileName({'Channel'},{curChannel},'fishImg','all');
    for j=1:numel(curFileList) 
        disp('Processing file:')
        disp(curFileList{j})
             
        %folder containing cropped loci & loc3
        [~,f,~] = fileparts(curFileList{j});
        cropFolder2 = fullfile(params.InputFolders.CropFolder,[f,'_crop']);
        % find the conditions in the current file
        fl = fs.getFileName({},{},'fishImg','all');
        idx = find( ismember(fl,curFileList{j}) );
        curCond = {};
        for k = 1:numel(fs.conditions)
            curCond = [curCond,fs.fList.(fs.conditions{k})(idx)];
        end
        
        if params.MatchLociSettings.testing
            ReferChannel = fs.getFileName({'Channel', 'FOV'},...
                {params.channelDescription.ReferChannel, cell2mat(curCond(2))},'fishImg');
        else
            ReferChannel = fs.getFileName({'Channel', 'Condition', 'FOV'},...
                {params.channelDescription.ReferChannel, curCond(2), cell2mat(curCond(3))},...
                'fishImg');
        end
        
        [~,fr,~] = fileparts(cell2mat(ReferChannel));
        cropFolder1 = fullfile(params.InputFolders.CropFolder,[fr,'_crop']);
        
        %load loc3 file
        loc1 = load(fullfile(cropFolder1,'location_mask_sorted.loc3'),'-ascii');
        loc2 = load(fullfile(cropFolder2,'location.loc3'),'-ascii');
        
        if isempty(loc1) || isempty(loc2)
            disp('Loc3 file is empty!')
            match = [];
            matchFileName = fullfile(params.Output.outFolder,'loci_match',...
                strcat('C',string(params.channelDescription.ReferChannel),'-',string(curChannel),fr(3:length(fr)),'.dloc3'));
            save(matchFileName,'match','-ascii');
            continue
        end
        
        %load cell mask
        disp(['Assigning cell ID for ', f]);
        MaskFile = fullfile(params.InputFolders.CellMaskFolder,strcat('C', string(DAPI),f(3:length(f)),'max_Mask.tiff'));
        Mask = readTifStackWithImRead(MaskFile);
        
        %sort loci into cell masks
        loc2 = sort_loci_into_cell_masks(loc2,Mask);
        save(fullfile(cropFolder2,'location_mask_sorted.loc3'),'loc2','-ascii');
        
        % convert to nm (assumes 3D)
        loc3nm1 = convert_loc_pix_to_nm(loc1,voxSize);
        loc3nm2 = convert_loc_pix_to_nm(loc2,voxSize);
        
        % optional: calculate dis hist between C1 and C2
        if compute_hist
            
            % compute distance matrix
            [dr, ~] = compute_distMatrix_from_loc(loc3nm1,loc3nm2,nDims);
        
            % compute pair distance histogram
            [n,x] = hist(dr( ~isnan(dr) ), randBinsCenters);
        
            % normalize
            fh = figure('name','Normalized pair distances between spots');
            plot(x,n./nRand);
            xlabel('Distance in nm');
            ylabel('Number of spot pairs');
        
            % save histogram
            if ~ exist(fullfile(params.Output.outFolder,'PairDistHist'),'dir')
                mkdir(fullfile(params.Output.outFolder,'PairDistHist'));
            end
            histFileName = fullfile(params.Output.outFolder,'PairDistHist',...
            strcat('C',string(params.channelDescription.ReferChannel),'-',string(curChannel),fr(3:length(fr)),'_pairDisHist.fig'));
            savefig(fh,histFileName);
            close(fh);
        end
              
        %loop through loci in C1
        missed1 = 0;
        match = [];
        for k = 1:length(loc1(:,1))
            cellID = round(loc1(k,7),12);
            %find loci in C2 that's in the same cell
            [idx,~] = find(round(loc2(:,7),12)==cellID);
            if isempty(idx)
                missed1 = missed1 + 1;
                fprintf('No match found in channel %d for loci %d in cell %d in channel %d.\n', ...
                    curChannel,k,cellID, params.channelDescription.ReferChannel);
            else
                %loop through loci in C2 and calculate dis
                dis = [];
                for l = 1:length(idx)
                    dis = [dis,...
                        sqrt((loc3nm1(k,1)-loc3nm2(idx(l),1)).^2+...
                        (loc3nm1(k,2)-loc3nm2(idx(l),2)).^2+...
                        (loc3nm1(k,3)-loc3nm2(idx(l),3)).^2)];
                end
                [dist, indx] = min(dis);
                if dist > params.MatchLociSettings.distThreshold
                    fprintf('No match found in channel %d for loci %d in cell %d in channel %d.\n', ...
                        curChannel,k,cellID, params.channelDescription.ReferChannel);
                    missed1 = missed1 + 1;
                else
                    delta = loc3nm1(k,1:3)-loc3nm2(idx(indx),1:3);
                    match = [match;[k,idx(indx),convert_loc_pix_to_nm(ceil(loc1(k,1:3))-ceil(loc2(idx(indx),1:3)),voxSize),...
                        dist,cellID,delta]];
                end
            end
        end
        if isempty(match)
            continue
        end
        %compute number of missed spots in C2
        missed2 = 0;
        match2=[];
        for k = 1:length(loc2(:,1))
            [idx,~] = find(match(:,2)==k);
            if length(idx) > 2
                fprintf('Loci %d in %s was matched for %d times.\n',k,f,length(idx));
                disp('Only the nearest two were kept')
                indx=[];
                for l = [1,2]
                    [~,indx1] = min(match(idx,4));
                    indx = [indx,indx1];
                end
                missed1 = missed1 + length(idx)-2;
                match2=[match2;match(indx,:)];
            else
                if isempty(idx)
                    missed2 = missed2+1;
                else
                    match2=[match2;match(idx,:)];
                end                
            end
        end
        match = match2;
        
        fprintf('Matched %d loci between %s and %s.\n', length(match(:,1)),fr,f);
        fprintf('Missed %d in %s\n', missed1, fr); 
        fprintf('Missed %d in %s\n', missed2, f);
        
        % Save match result and add to file list
        matchFileName = fullfile(params.Output.outFolder,'loci_match',...
            strcat('C',string(params.channelDescription.ReferChannel),'-',string(curChannel),fr(3:length(fr)),'.dloc3'));
        fs.setFileName(fs.conditions,curCond,'loci_match',matchFileName);
        save(matchFileName,'match','-ascii');
        
        %generate a plot for visual check (plot on x y plane)
        %matched C1 loci
        x = loc1(match(:,1),1:2);
        scatter(x(:,1),x(:,2),20,'blue','filled');
        hold on
        %missed C1 loci
        idx = setdiff(1:length(loc1(:,1)),match(:,1));
        x = loc1(idx,1:2);
        scatter(x(:,1),x(:,2),20,'cyan','filled')
        hold on
        %match C2 loci
        x = loc2(match(:,2),1:2);
        scatter(x(:,1),x(:,2),20,'red','filled');
        hold on
        %missed C2 loci
        idx = setdiff(1:length(loc2(:,1)),match(:,2));
        x = loc2(idx,1:2);
        scatter(x(:,1),x(:,2),20,'magenta','filled')
        legend('C1 matched loci','C1 missed loci','C2 matched loci','C2 missed loci')
        hold off
        savefig(fullfile(params.Output.outFolder,'loci_match',strcat('C',string(params.channelDescription.ReferChannel),'-',string(curChannel),fr(3:length(fr)),'.fig')))
        close();
    end
end