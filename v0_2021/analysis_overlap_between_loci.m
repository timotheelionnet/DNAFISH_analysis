addpath('../fileSet');
addpath('../iniconfig');
addpath('../cbrewer');

%% get list of files
try
    configFileName;
catch
    % config file path
    configFileName = '20210203_DNA_FISH_overlap_config.ini';
end

% create file list object
fs = fileSet;
fs.buildFileSetFromConfig(configFileName);

% reformat channel, replicate and FOV as numbers
convertNonNumeralStrings = 1; 
fs.toNumeric({'Channel','TechnicalReplicate','FOV'},convertNonNumeralStrings);

%read parameters from config file relative to analysis into a parameters object
params = readParamsFromIniFile(configFileName);

%% Calculate overlap between pairs
voxsize = [params.Settings.voxSize_dxy,...
    params.Settings.voxSize_dxy,...
    params.Settings.voxSize_dz];

for i = 1:numel(fs.fList.(width(fs.fList)))
    %load match file
    match = load(fs.fList.(width(fs.fList)){i},'-ascii');
    if isempty(match)
        continue
    end
    [~,f,~] = fileparts(fs.fList.(width(fs.fList)){i});
    disp('Processing file:')
    disp(f)
    C1 = fs.fList.(1){i};
    C2 = fs.fList.(2){i};
    CropFolder1 = fullfile(params.InputFolders.CropFolder,...
            strcat('C',C1,f(5:length(f)),'_crop'));
    CropFolder2 = fullfile(params.InputFolders.CropFolder,...
            strcat('C',C2,f(5:length(f)),'_crop'));
    %find out nDigitsMax for each image, test up to 5 digits
    for l = 1:5
        try 
           locus1 = timtiffread(fullfile(CropFolder1,...
                ['locus_bgcorr',add_digits(1,l),'.tif']));
        catch
           continue
        end
        nDigitsMax1 = l;
        break
    end
    
    for l = 1:5
        try 
           locus2 = timtiffread(fullfile(CropFolder2,...
                ['locus_bgcorr',add_digits(1,l),'.tif']));
        catch
           continue
        end
        nDigitsMax2 = l;
        break
    end
        
    overlap = [];
    for j = 1:length(match(:,1))
        fprintf('Processing %s C%s locus%d and C%s locus%d...\n',...
            f(6:length(f)),C1,match(j,1),C2,match(j,2))
         
        %load cropped images
        locus1 = readTifStackWithImRead(fullfile(CropFolder1,...
            ['locus_bgcorr',add_digits(match(j,1),nDigitsMax1),'.tif']));
        locus2 = readTifStackWithImRead(fullfile(CropFolder2,...
            ['locus_bgcorr',add_digits(match(j,2),nDigitsMax2),'.tif']));
        %normalize by sum
        locus1n = double(locus1)/sqrt(sum(double(locus1).^2,'all'));
        locus2n = double(locus2)/sqrt(sum(double(locus2).^2,'all'));
        %initialize a big 'image' consist of zeros
        sizeb = size(locus1);  
        deltanm = match(j,3:5);
        delta = convert_loc_pix_to_nm(deltanm,1./voxsize);
        I = zeros(sizeb*2+abs(delta)*2);
        %put locus1 and locus2 in the image        
        loc1 = 1.+delta;
        loc2 = [1,1,1];%relative location of locus1 and 2
        for k = 1:3
            if delta(k)<1
                loc1(k) = 1;
                loc2(k) = 1-delta(k);%to make sure idx is positive
            end
        end
        I(loc1(1):loc1(1)+sizeb(1)-1,loc1(2):loc1(2)+sizeb(2)-1,loc1(3):loc1(3)+sizeb(3)-1) = locus1n;
        I(loc2(1):loc2(1)+sizeb(1)-1,loc2(2):loc2(2)+sizeb(2)-1,loc2(3):loc2(3)+sizeb(3)-1) =...
            I(loc2(1):loc2(1)+sizeb(1)-1,loc2(2):loc2(2)+sizeb(2)-1,loc2(3):loc2(3)+sizeb(3)-1)+locus2n;
        %calculate correlation
        corr = (sum(I.*I,'all')-sum(locus1n.*locus1n,'all')-sum(locus2n.*locus2n,'all'))/2;
        %save correlation
        overlap = [overlap; corr];
    end
    
    %save corr in file
    if ~ exist(params.Output.outFolder,'dir')
        mkdir(params.Output.outFolder);
    end
    save(fullfile(params.Output.outFolder,strcat(f,'_overlap.txt')),'overlap','-ascii')
end


%%
function s = add_digits(n,nDigitsMax)
    nDigits = ceil(log(n+1)/log(10));
    nZerosToAdd = nDigitsMax - nDigits;
    s = [repmat('0',1,nZerosToAdd),num2str(n)];
end