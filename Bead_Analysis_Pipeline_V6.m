%% Bead Analysis Pipeline V6
% Updated 12/5/23 to read loc4 files. 
% Lionnet Lab | Finn Clark
% For modeling channel dependent displacement in multi-ch imaging systems
% Input = path to bead scan filepath
% Output = NN analysis and displacement 
% field fit functions



%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Define User Settings~~~~~~~~~~~~~~~~~~~~~~~~

% ~~~~~~You must fill this out~~~~~~~~~~~~~~~~~~

%Make a fileset object to store filepaths in a table
fs = fileSet;

% beadConfigPath = "H:\Finn_DriveH\20240226_PPARG_Staggered_Timecourse\beads\session3_beads.ini";
% Do you want to run AIRLOCALIZE: 1
% loc4 files already saved in directory: 0
runAirlocalize = 0;
trainAirlocalize = 0;

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
    myFileLists(i) = {fs.getFileName({'Channel'},{allBeadChannels(i)},'fishImg','all')};   

    disp(append("You have loaded ", string(height(myFileLists{i})), ...
    " Images for Channel " ,string(allBeadChannels(i))))
end
% Print ref channel
refChannel= myParams.channelDescription.refChannel;
disp(append("Channel ", string(refChannel), " is your reference channel for registration"));




%% Run AIRLOCALIZE if necessary


if runAirlocalize ==1
    if trainAirlocalize ==1
    
        % -----------------------------GUI for getting PSF and thresh
    
        
        %I want to go through and manually get fit params using thefirst image of
        %each channel.
         
        locpar = cell(1, numel(allBeadChannels));
        for i = 1:numel(allBeadChannels)
            disp('Processing');
            disp(myFileLists{i}{1});
            [~, par4_path] = perform_detection_on_single_img_once_Finn(myFileLists{i}{1}, 0, 1, 1,2); %Manual training on first image for ith channel
            
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
%      for i = 1:1
        for j=1:numel(myFileLists{i})
            
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
            fl = fs.getFileName({},{},'fishImg','all');
            idx = find(ismember(fl,myFileLists{i}{j}));
            
          
            % curCond = {};
            % 
            % for k = 1:numel(fs.conditions)
            %     curCond = [curCond,fs.fList.(fs.conditions{k})(idx)];
            % end
            % %Save the images loc4 filename in the fileSet
            % fs.setFileName(fs.conditions,curCond,'loc4',locFileName);
            
            %Now save the actual file
            save(locFileName,'loc','-ascii');
            
    
        end
    end
disp('Airlocalize complete');
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


myParams = readParamsFromIniFile(beadConfigPath);
%load any FOV to get image size
curChannel = myParams.channelDescription.fishChannelList(1);
curImgFileName = fs.getFileName({'Channel'},{curChannel},'fishImg','first');
curImgFileName = curImgFileName{1}; %extract filepath as string from cell
imSize = getImageSizeFromFile(curImgFileName); % get image size in pixels

% create a directory for results output 
resDir = fullfile(fileparts(fileparts(curImgFileName)), 'res');


voxSize = [myParams.Settings.voxSize_dxy,...
            myParams.Settings.voxSize_dxy,...
            myParams.Settings.voxSize_dz];


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
    
    
    for j = 1:numel(myRegFl) 
        %load jth FOV and convert to nm REMEMBER FOVS are in order of
        %leading digit low to high!!
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
curImgFileName = fs.getFileName({'Channel'},{curChannel},'fishImg','first');
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
curImgFileName = fs.getFileName({'Channel'},{curChannel},'fishImg','first');
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
xlim([0,1000]);

savefig(resid_cdf, fullfile(resDir, 'cdfs'));

%%

%Manually plot and fit a camera and dimension
myDim=1;
myCam = 1; 

%reg chan 1, x,y, dx
myX = allNn_data{myCam}(:,4);
myY = allNn_data{myCam}(:,5);
myZ = allNn_data{myCam}(:,7+myDim); %dr for given dim

myZ1 = allNn_data{myCam}(:,8);
myZ2 = allNn_data{myCam}(:,9);

[myfit, gof, output] = fit([myX, myY], myZ, fitType, 'Normalize', 'on');

myDims = ["x" "y" "z"];
myCams = ["1" "2"];

%PLOT THE FIT dr as a function of reg x,y
figure;
plot(myfit, [myX, myY], myZ);
title(append("d" , myDims(myDim), " for camera ", myCams(myCam)));
xlabel('X (nm)');
ylabel('Y (nm)');
zlabel('di (nm)');

%get predicted displacement
zFitted = feval(myfit, myX, myY);
%Get displacement residual s
zResids = zFitted - myZ; 


%plot residual map
zResidsMap = min(quantile(zResids,0.9), max(quantile(zResids,0.1),zResids));
figure;
scatter3(myX, myY, zResids, 12*ones(size(myX)),zResidsMap,'filled');
% quiver(myX, myY, myZ1, myZ2); %plots vectors originating at myX myY and with direction myZ1 myZ2
title(append("XY displacement vec field for camera ", myCams(myCam)));
xlabel('X (nm)');
ylabel('Y (nm)');

%%




disp('Bead Analysis Complete');