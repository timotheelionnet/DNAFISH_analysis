%% Bead Analysis Pipeline V6 MODULAR
% Updated 3/14/24
% Lionnet Lab | Finn Clark
% For modeling channel dependent displacement in multi-ch imaging systems
% Input = path to bead scan filepath
% Output = NN analysis and displacement 
% field fit functions

% clear
%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define User Settings~~~~~~~~~~~~~~~~~~~~~~~~

%lets check for an ini file

% ~~~~~~You must fill this out~~~~~~~~~~~~~~~~~~

%Make a fileset object to store filepaths in a table
fs = fileSet;

% beadConfigPath = "H:\Finn_DriveH\20240226_PPARG_Staggered_Timecourse\beads\session2_beads.ini";
beadConfigPath = cfg_path;

% Do you want to run AIRLOCALIZE: 1
% loc4 files already saved in directory: 0
% runAirlocalize = 1;
trainAirlocalize = 1;

% Do you want to define an ROI within which to construct your model?
% I normally do not use this feature. We can always set an ROI when we
% apply the model to real data. That way edge effects of fitting are kept
% at the very edges of the FOV where little data will be present

crop = 0; % 1 Yes or 0 No
cropRad = 0.5; % Only if crop =1. ROI will be a square centered on the FOV with sidelength = 2* cropRad * (sidelength of image);
zCropStDevs = 6; % How many stdevs in z from the mean do we allow a point to be. Normally doesn't matter. 

% Set the distance threshold for matching nearest neighbors between
% channels in nm
myNnThresh = 1000;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




% convert to char for downstream compat
beadConfigPath = char(beadConfigPath);
% build fs from the beadscan cfg
fs.buildFileSetFromConfig(beadConfigPath);


% reformat channel, replicate and FOV as numbers
%This must be done in order to use the getFileName function!!!!
convertNonNumeralStrings = 1; 
fs.toNumeric({'Channel','FOV'},convertNonNumeralStrings);

% Report fList size
numF = numel(unique(fs.fList.FOV));
numCh = numel(unique(fs.fList.Channel));
disp("~~~~~~Loaded Fileset Containing~~~~~~");
disp(string(numF) +" FOVs");
disp(string(numCh) + " Channels");
disp(" ");

%% Params from fileset

%Add a column to your fileSet object to hold localization data

%Read the parameters for your imaging session/data
myParams = readParamsFromIniFile(beadConfigPath);

%Extract ist of your channels with data to be localized
allBeadChannels = myParams.channelDescription.fishChannelList;

%Initialize cell array that will contain a column of image filepaths for each
%channel. 
myFileLists = cell(1, numel(allBeadChannels)); 

for i=1:numel(allBeadChannels) %index based on how many fish channels
    myFileLists(i) = {fs.getFileName({'Channel'},{allBeadChannels(i)},'beadImg','all')};   

    disp(append("You have loaded ", string(height(myFileLists{i})), ...
    " Images for Channel " ,string(allBeadChannels(i))))
end
% Print ref channel
refChannel= myParams.channelDescription.refChannel;
disp(append("Channel ", string(refChannel), " is your reference channel for registration"));



% create a directory for results output 
resDir = fullfile(bead_dir, 'res');

if ~exist(resDir, 'dir')
    mkdir(resDir)
    disp('Bead analysisi results will be saved to')
    disp(resDir)
end

%% Run AIRLOCALIZE if necessary



if runAirlocalize == 1
    % collect all loc in a cell chwere cell n is  bead channel n
    compiled_loc = [];
    if trainAirlocalize ==1
    
        % -----------------------------GUI for getting PSF and thresh
    
        
        %I want to go through and manually get fit params using thefirst image of
        %each channel.
         
        locpar = cell(1, numel(allBeadChannels));
        for i = 1:numel(allBeadChannels)
        % for i = 1:1 

            curImPath = myFileLists{i}{1}; % get first in list

            disp('Processing');
            disp(curImPath);

            % get approx of how many beads are in the FOV

            n_beads = count_binary_beads_for_pipe(curImPath);

            % try running the optimizer
            cur_thresh = 1; % initial guess for SD thresh val

            bead_AL_thresh_optimizer % run
            
            % input = img path, interactive psf, interactive thresh,
            [~, par4_path] = perform_detection_on_single_img_once_Finn(myFileLists{i}{1}, 0, 0, 1,2, 'SD', cur_thresh); %Manual training on first image for ith channel
            
            locpar{i} = par4_path; %Store first par4 file for ith channel
           
            
        end
    end

    %----------------------------Par loop to process remaining FOVS

    fs.addFileType('loc',[]);
    %Now I have trained AIRLOCALIZE on all channels and stored the fit params
    %Time to loop through all of my images and analyze those using previosuly
    %stored fit params

    imageCount= append("ch:", string(numel(allBeadChannels)), " file:", string(numel(myFileLists{1})));
    disp("hi");
    
    for i = 1:numel(allBeadChannels)
    % for i = 1:1
        for j=1:numel(myFileLists{i})

            cur_ch = allBeadChannels(i);
            
            disp(append("processing file ch:", string(i), " file:", string(j), " out of ", imageCount));
    
            
            %Using ith channel parameters, replace data with jth image in
            %channel i
            % locpar{i}.data = readTifStackWithImRead(myFileLists{i}{j});
    
            % run airlocalize without GUI for the rest of the files
            loc = parlooper(myFileLists{i}{j}, locpar{i});
            
            disp(append("Processed", string(myFileLists{i}{j})));
            
            %For each image we want to store loc4 data and its file name in the
            %correct row of our loc4 fileSet column. loc4 is overwritten each
            %iteration
            % loc4 = res{i}.final_pix;
    
            % build the loc filename to enter in the fileSet list
            [d,f,~] = fileparts(myFileLists{i}{j});
            locFileName = fullfile(d,[f,'.loc4']);
            
            % find the conditions of the current file
            % fl = fs.getFileName({},{},'fishImg','all');
            % idx = find(ismember(fl,myFileLists{i}{j}));
            
          
            % curCond = {};
            % 
            % for k = 1:numel(fs.conditions)
            %     curCond = [curCond,fs.fList.(fs.conditions{k})(idx)];
            % end
            % %Save the images loc4 filename in the fileSet
            % fs.setFileName(fs.conditions,curCond,'loc4',locFileName);

            % compile in cell array
            compiled_loc = vertcat(compiled_loc, [repmat(allBeadChannels(i), [height(loc),1]), loc(:,1:3), repmat(j, [height(loc),1]) ] );
            
            %Now save the actual file
            save(locFileName,'loc','-ascii');
            
    
        end
    end
disp('Airlocalize complete');

% convert compiled loc to table format
% saves in units of pix
compiled_loc_tab = array2table(compiled_loc);
compiled_loc_tab.Properties.VariableNames = {'channel', 'x', 'y', 'z', 'FOV'};

save_name = fullfile(resDir, 'compiled_loc_tab.csv');
writetable(compiled_loc_tab, save_name);

disp('Compiled loc table saved to:')
disp(save_name)
else
    disp('Not running Airlocalize');
end




%% Reload file list object to make sure that locs are loaded
fs = fileSet;
configFileName = beadConfigPath;
fs.buildFileSetFromConfig(configFileName);


% reformat channel, replicate and FOV as numbers
convertNonNumeralStrings = 1; 
fs.toNumeric({'Channel','FOV'},convertNonNumeralStrings);

%% All by all dist hist between non ref and ref channel

% compile loc if not running AL (already done)
if runAirlocalize == 0
    compiled_loc = []; % channel, x,y,z, FOV index (j)
end

myParams = readParamsFromIniFile(beadConfigPath);
%load any FOV to get image size
curChannel = allBeadChannels(1);
curImgFileName = fs.getFileName({'Channel'},{curChannel},'beadImg','first');
curImgFileName = curImgFileName{1}; %extract filepath as string from cell
imSize = getImageSizeFromFile(curImgFileName); % get image size in pixels



% voxSize = [myParams.Settings.voxSize_dxy,...
%             myParams.Settings.voxSize_dxy,...
%             myParams.Settings.voxSize_dz];


%Dimensionality of your images
nDims = numel(voxSize);

imSize_nm = convert_loc_pix_to_nm(imSize,voxSize);

% this bin size always works so hard coding it in
distMat_hist_binSize = 50;
nRandSpots = 10000; % number of spots used for pair dist normalization



%Generate random spots for normalizing pair distance in bead images 
[nRand,randBinsCenters] = generate_pairDist_hist_randSpots(...
    nRandSpots,...
    imSize_nm,...
    distMat_hist_binSize);



% Now I want to get a file list for each channel (Taken from chunk above).

% Define vars for refrence channel and non reference channels
refCh = myParams.channelDescription.refChannel;
regChannels = setdiff(allBeadChannels, myParams.channelDescription.refChannel);

%Get cell array containing all non reference channel file names
regChLocLists = cell(1, numel(regChannels));
for i = 1:numel(regChannels)
    
    regChLocLists{i} = fs.getFileName({'Channel'}, {regChannels(i)} ,'loc','all');
     
end    

% load all my reference channel loc4 paths in a fileList

reflocList = fs.getFileName({'Channel'},{refCh},'loc','all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Working dist hist FOV by FOV

pairHistCounts = cell(1, numel(regChannels)); %Preallocate one cell per registration channel
for i = 1:numel(regChannels)
    %load ith non reference channel as a fileList
    myRegFl = regChLocLists{i};

    cur_reg_ch = regChannels(i);
    
    
    for j = 1:numel(myRegFl) 
        %load jth FOV and convert to nm REMEMBER FOVS are in order of
        %leading digit low to high!!
        locRef = load(reflocList{j});
        locReg = load(myRegFl{j});

        % compile loc if not running AL
        if runAirlocalize == 0
        
            % take the ref loc the first time
            if i ==1
                % channel, x,y,z, FOV index (j)
                compiled_loc = vertcat( compiled_loc, [repmat(refCh, [height(locRef),1]), locRef(:, 1:3),  repmat(j, [height(locRef),1]) ]);
            end
            % channel, x,y,z, FOV index (j)
            compiled_loc = vertcat( compiled_loc, [repmat(cur_reg_ch, [height(locReg),1]), locReg(:, 1:3),  repmat(j, [height(locReg),1]) ]);
    

        end

        locnmRef = convert_loc_pix_to_nm(load(reflocList{j}), voxSize);
        locnmReg = convert_loc_pix_to_nm(load(myRegFl{j}), voxSize);
        
              
        %compute distance matrix for jth FOV between ith reg channel and ref
        % channel
        [dr, ~] = compute_distMatrix_from_loc(locnmRef,locnmReg, nDims);
        
        % compute pair distance histogram. n is counts ber bin
        [n,x] = hist(dr( ~isnan(dr) ), randBinsCenters);
        
        %store the normalized counts for the ith channel and jth FOV
        pairHistCounts{i} = vertcat(pairHistCounts{i}, n ./ nRand);
        
        
        
        
    end
end   

if runAirlocalize == 0
    % save compiled loc (pixel units)
    compiled_loc_tab = array2table(compiled_loc);
    compiled_loc_tab.Properties.VariableNames = {'channel', 'x', 'y', 'z', 'FOV'};

    save_name = fullfile(resDir, 'compiled_loc_tab.csv');
    writetable(compiled_loc_tab, save_name);
   
    disp('Compiled loc table saved to:')
    disp(save_name)
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot and save all by all distance plot
histDir = fullfile(resDir, 'DistHists');
if ~ exist(fullfile(histDir),'dir')
    mkdir(fullfile(histDir));
end

%%Working plot function using boolean mask to exclude nans and infs
idx = ~isnan(pairHistCounts{1}(1,:)) & ~isinf(pairHistCounts{1}(1,:)); %gets problematic column indices for 1st ch
for i=1:size(pairHistCounts,2) % loop through channels 
    f = figure;
    plot(randBinsCenters(idx),mean(pairHistCounts{i}(:,idx))); %plto avg bin val across all FOVs
    xlabel("Pair Distances (nm)");
    ylabel("Counts");
    title(append("Pair Dist Hist: Ch", string(regChannels(i)), " vs Reference"));
    xlim([0 0.5e4])

    saveName = fullfile(histDir, 'PairDistHist_Ch' + string(regChannels(i)) + ' vs ref');
    savefig(f ,saveName);
    close(f);
    
end
disp("Saved pair dist hists to: " + string(histDir));


%% Construct NN matrix

refCh = myParams.channelDescription.refChannel;
reflocList = fs.getFileName({'Channel'},{refCh},'loc','all');
regChLocLists = cell(1, numel(regChannels));

for i = 1:numel(regChannels)
    
    regChLocLists{i} = fs.getFileName({'Channel'}, {regChannels(i)} ,'loc','all');        
end    

% Using the get nn matrix func
% inputs are ref loclist, cell array of reg loc lists, nn thresh in
% relevant units, nDims, and vox size in desired units per pixel in xyz
% Output is [refx, refy, refx], [regx, regy, regz], 3d dist, xdist,ydist, zdist
% as well as percentage of spots that did not get matched to an NN. Often
% due to different thresholdng in AIRLOCALIZE

allNn_data = get_nn_matrix3(reflocList, regChLocLists, myNnThresh, nDims=3, voxSize=voxSize);



%% Crop FOV if desired

if crop == 1

    %define center of ROI
    center = [imSize_nm(1)/2 , imSize_nm(2)/2];
    %radius from center
    rad = center(1)*cropRad;
    roiIDX = cell(1, numel(allNn_data));
    
    %loop through each registration channel and crop by channel
    for i =1:numel(allNn_data)
        
        myNN_mat = allNn_data{i};

        
        %get mask for where x and y and Z are within the ROI
        roiIDX{i} = myNN_mat(:,1) > (center(1) - rad)  & myNN_mat(:,1) < (center(1) + rad) ...
            & myNN_mat(:,2) > (center(2) - rad) & myNN_mat(:,2) < (center(2) + rad) ...
            & myNN_mat(:,3) > (  mean(myNN_mat(:,3)) - 6*std(myNN_mat(:,3))  )...
            & myNN_mat(:,3) < (  mean(myNN_mat(:,3)) + 6*std(myNN_mat(:,3))  );
        
        
        allNn_dataROI{i} = allNn_data{i}(roiIDX{i}, :);
    end
    
    disp('ROI Applied');

else
    disp('No Crop Applied')

end

%% Plot Scatter of our bead localizations in 3 channels
disp('Plotting bead loc data');

% get an example file name 
curChannel = myParams.channelDescription.fishChannelList(1);
curImgFileName = fs.getFileName({'Channel'},{curChannel},'beadImg','first');
curImgFileName = curImgFileName{1}; %extract filepath as string from cell

% If we are cropping, we must filter the NN data 
if  crop ==1 
    allNn_data = allNn_dataROI;
end

% plot all bead data
x_ref = [];
fScat = figure ; 
legTxt = ['Ch'+ string(refCh) + ' (ref)'];

for i = 1:numel(allNn_data)
    x_ref = allNn_data{i}(:,1:3);
   
    scatter3(x_ref(:,1),x_ref(:,2), x_ref(:,3), '.');
    
    hold on
end

for i = 1:numel(allNn_data)
    
    x_reg = allNn_data{i}(:,4:6);
    
    
    legTxt = [legTxt, 'Ch'+ string(regChannels(i))];
    
    scatter3(x_reg(:,1),x_reg(:,2), x_reg(:,3), '.','DisplayName', 'Ch'+ string(regChannels(i)));
    hold on
    
end


title('Matched Nearest Neighbors');
xlabel('X (nm)');
ylabel('Y (nm)');
zlabel('Z (nm)');
legend(legTxt);

% Save a scatter plot of the data
savefig(fScat, fullfile(resDir, 'Matched Nearest Neighbors Scatter Plot'));
close(fScat);

disp("Plotting complete");
disp("Saved scatter plot to: " + string(resDir));
%% Perform a polynomial fit to the displacement field and quantify accuracy


disp('Modelng displacement field');

% get an example file name for modularity
curChannel = myParams.channelDescription.fishChannelList(1);
curImgFileName = fs.getFileName({'Channel'},{curChannel},'beadImg','first');
curImgFileName = curImgFileName{1}; %extract filepath as string from cell

% choosing a fit type. polyxy means polynomial fit function degree x in x
% and degree y in y
fitType = 'poly33';
myFits= cell(max(allBeadChannels));

%going to compile dr before and after correction
preReg_dr= cell(size(regChannels));
postReg_dr= cell(size(regChannels));

disp('Fit type: '+ string(fitType));
for i = 1:numel(regChannels)
    
    r_reg = [allNn_data{i}(:,4) allNn_data{i}(:,5) allNn_data{i}(:,6)]; 
    r_ref = [allNn_data{i}(:,1) allNn_data{i}(:,2) allNn_data{i}(:,3)];
    
    dr = r_reg - r_ref; 
         
    p = [r_reg(:, 1) r_reg(:, 2)]; %These are my XY coords in the reg channel
    v = [dr(:,1) dr(:,2) dr(:,3)]; %Store dx dy dz
    
    %Generate a poly fit for dx dy dz as a function of XY
    [dxFit, dxgof, dxoutput] = fit(p, v(:,1),fitType, 'Normalize', 'on');
    [dyFit, dygof, dyoutput] = fit(p, v(:,2),fitType, 'Normalize', 'on');
    [dzFit, dzgof, dzoutput] = fit(p, v(:,3),fitType, 'Normalize', 'on');

    %Store fit objects, ith cam, displacement in jth dimension {fit object,
    %gof params, output params}
    myFits{regChannels(i), refCh}{1} = {dxFit, dxgof, dxoutput};
    myFits{regChannels(i), refCh}{2} = {dyFit, dygof, dyoutput};
    myFits{regChannels(i), refCh}{3} = {dzFit, dzgof, dzoutput};
    
    
    %Evaluate fit function at all XY coords for dx dy dz
    fEvals{i}{1} = feval(dxFit, p);
    fEvals{i}{2} = feval(dyFit, p);
    fEvals{i}{3} = feval(dzFit, p);
    
    %Get residuals of form empirical - fitted
    %reg channel i and dx dy dz 
    fResids{i}{1} = v(:,1) - fEvals{i}{1};
    fResids{i}{2} = v(:,2) - fEvals{i}{2};
    fResids{i}{3} = v(:,3) - fEvals{i}{3};
    
%     remove BIG RESIDS

    % fResids{i}{1} = fResids{i}{1}(abs(fResids{i}{1}) < quantile(abs(fResids{i}{1}), 0.95));
    % fResids{i}{2} = fResids{i}{2}(abs(fResids{i}{2}) < quantile(abs(fResids{i}{2}), 0.95));
    % fResids{i}{3} = fResids{i}{3}(abs(fResids{i}{3}) < quantile(abs(fResids{i}{3}), 0.95));
  

    % store nn disp before and after correction
    preReg_dr{i} = v;
    postReg_dr{i} = [fResids{i}{1}, fResids{i}{2}, fResids{i}{3}];

    % Get summary stats for each set of resids 
    % cell i = reg channel i 
    % each cell contains cols: rmsErr, medAbsErr, sdErr for rows: dx dy dz
    
    for j = 1:nDims
        [rmsqErr, medAbsErr, sdErr] = deal(rms(fResids{i}{j}), median(fResids{i}{j}), std(fResids{i}{j}));
        
        summary{1,i}(j+1,:) = [string(rmsqErr), string(medAbsErr), string(sdErr)];
    end
    summary{1,i}(1,:) = ["rmsqErr", "medAbsErr", "sdErr"];
   
    
%     [v_unique, p_unique] = groupsummary(v,p, @max);
    
      
end

% Display summaries for arbitrary number of channel registrations
for i = 1:numel(regChannels)
    myChan = regChannels(i);
    disp("~~~~Channel " + string(myChan) +" Reg Model Summary~~~~");
    disp(summary{i});
end

% Save the fits as a matlab object
[~,samp_name] = fileparts(fileparts(resDir));
save(fullfile(resDir, strcat(samp_name, '_fits')),'myFits');
disp("Saved registration models in " + string(resDir));

%% Visualizing

% chNames = myParams.channelDescription.fishChannelDescription;
% lets visualize registration correction 
resid_cdf = figure;

res_legTxt = [];
% go thru all i registered channels
for i = 1:numel(fEvals)
    
    % pre reg data
    dv = preReg_dr{i};
    preReg_dr_mag = sqrt(dv(:,1).^2 + dv(:,2).^2 + dv(:,3).^2) ;

    cdfplot(preReg_dr_mag);
    hold on

    res_legTxt = [res_legTxt, 'Raw Ch' + string(regChannels(i))];
    
    % post reg data
    dv = postReg_dr{i};
    postReg_dr_mag = sqrt(dv(:,1).^2 + dv(:,2).^2 + dv(:,3).^2) ;

    cdfplot(postReg_dr_mag);
    hold on

    res_legTxt = [res_legTxt, 'Reg Ch' + string(regChannels(i))];
end
legend(res_legTxt);

title('Registration Error for Bead Corrected Channels');
xlabel('3D Displacement (nm)')
% xlim([0,1000]);

savefig(resid_cdf, fullfile(resDir, 'goodnessOfFit_cdfs'));

%%



disp('Bead Analysis Complete');