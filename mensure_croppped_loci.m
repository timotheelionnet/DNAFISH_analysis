addpath('../fileSet');
addpath('../iniconfig');

%% get list of files

% config file path
configFileName = '20210203_DNA_FISH_mensure_config.ini';

% create file list object
fs = fileSet;
fs.buildFileSetFromConfig(configFileName);

% reformat channel, replicate and FOV as numbers
convertNonNumeralStrings = 1; 
fs.toNumeric({'Channel','TechnicalReplicate','FOV'},convertNonNumeralStrings);

%read parameters from config file relative to analysis into a parameters object
params = readParamsFromIniFile(configFileName);

%% calculate centriods 

%loop through crop folders
for i=1:numel(params.channelDescription.fishChannelList)
        curChannel = params.channelDescription.fishChannelList(i);
        curFileList = fs.getFileName({'Channel'},{curChannel},'fishImg','all');
        for j=1:numel(curFileList)
            
            %load cropped images
            [~,f,~] = fileparts(curFileList{j});
            CropFolder = fullfile(params.InputFolders.CropFolder,strcat(f,'_crop'));
            disp('Processing cropped images of:')
            disp(f)
    
            %find out nDigitsMax for each image, test up to 5 digits
            for l = 1:5
                try 
                locus1 = timtiffread(fullfile(CropFolder,...
                     ['locus_bgcorr',add_digits(1,l),'.tif']));
                catch
                    continue
                end
                nDigitsMax1 = l;
                break
            end
    
            %load loc3 file
            loc = load(fullfile(CropFolder,'location.loc3'),'-ascii');
    
            %loop through loci and calculate centriods
            cen = [];
            for k = 1:length(loc(:,1))
                %load cropped image
                locus = readTifStackWithImRead(fullfile(CropFolder,...
                    ['locus_bgcorr',add_digits(k,nDigitsMax1),'.tif']));
        
                if i == 1 && j == 1 && k == 1
                    %generate meshgrid
                    size0 = size(locus);
                    [x,y,z] = meshgrid(1:size0(1),1:size0(2),1:size0(3));
                end
        
                %calculate centroid
                xc = sum(x.*locus,'all')/sum(locus,'all');
                yc = sum(y.*locus,'all')/sum(locus,'all');
                zc = sum(z.*locus,'all')/sum(locus,'all');
        
                cen = [cen;[xc,yc,zc]];       
            end
            %save centroid
            FileName = fullfile(CropFolder,'centroid.loc3');
            save(FileName,'cen','-ascii');  
        end
end

%% calculate radio of gyration of loci
%loop through crop folders
for i=1:numel(params.channelDescription.fishChannelList)
    curChannel = params.channelDescription.fishChannelList(i);
    curFileList = fs.getFileName({'Channel'},{curChannel},'fishImg','all');
    for j=1:numel(curFileList)
        
        %load cropped images
        [~,f,~] = fileparts(curFileList{j});
        CropFolder = fullfile(params.InputFolders.CropFolder,strcat(f,'_crop'));
        disp('Processing cropped images of:')
        disp(f)
    
        %find out nDigitsMax for each image, test up to 5 digits
        for l = 1:5
            try 
            locus1 = timtiffread(fullfile(CropFolder,...
                    ['locus_bgcorr',add_digits(1,l),'.tif']));
            catch
                continue
            end
            nDigitsMax1 = l;
            break
        end
    
        %load centroid file
        cen = load(fullfile(CropFolder,'centroid.loc3'),'-ascii');
    
        %loop through loci and calculate Rg
        r = [];
        for j = 1:length(cen(:,1))
        
            %load cropped image
            locus = readTifStackWithImRead(fullfile(CropFolder,...
                ['locus_bgcorr',add_digits(j,nDigitsMax1),'.tif']));
        
            %generate meshgrid if not already exist
            try
                [x,y,z];
            catch
                size0 = size(locus);
                [x,y,z] = meshgrid(1:size0(1),1:size0(2),1:size0(3));
            end
        
            %centroid coordinate
            cen0 = cen(j,:);
            xc = cen0(1);
            yc = cen0(2);
            zc = cen0(3);
            
            %calculate Rg
            rx = sum(sqrt(power((x-xc),2).*locus),'all')/sum(locus,'all');
            ry = sum(sqrt(power((y-yc),2).*locus),'all')/sum(locus,'all');
            rz = sum(sqrt(power((z-zc),2).*locus),'all')/sum(locus,'all');
        
            r = [r;[rx,ry,rz]];       
        end
        %save Rg
        FileName = fullfile(CropFolder,'Rg.loc3');
        save(FileName,'r','-ascii');  
    end
end

%% calculate volumn of loci by threshold

for i=1:numel(params.channelDescription.fishChannelList)
    curChannel = params.channelDescription.fishChannelList(i);
    curFileList = fs.getFileName({'Channel'},{curChannel},'fishImg','all');
    for j=1:numel(curFileList)

        %load cropped images
        [~,f,~] = fileparts(curFileList{j});
        CropFolder = fullfile(params.InputFolders.CropFolder,strcat(f,'_crop'));
        disp('Processing cropped images of:')
        disp(f)
    
        %find out nDigitsMax for each image, test up to 5 digits
        for l = 1:5
            try 
            locus1 = timtiffread(fullfile(CropFolder,...
                    ['locus_bgcorr',add_digits(1,l),'.tif']));
            catch
                continue
            end
            nDigitsMax1 = l;
            break
        end
    
        %load loc3 file
        loc = load(fullfile(CropFolder,'location.loc3'),'-ascii');
    
        %loop through loci and calculate volumn
        count = [];
        for j = 1:length(loc(:,1))
            %load cropped image
            locus = readTifStackWithImRead(fullfile(CropFolder,...
                ['locus_bgcorr',add_digits(j,nDigitsMax1),'.tif']));
        
            %count number of pixel above threshold
            idx = find(locus > params.Settings.threshold);
            count = [count;length(idx)];
        end
    
        %save volumn
        FileName = fullfile(CropFolder,'volumn.txt');
        save(FileName,'count','-ascii'); 
    end
end

%%
function s = add_digits(n,nDigitsMax)
    nDigits = ceil(log(n+1)/log(10));
    nZerosToAdd = nDigitsMax - nDigits;
    s = [repmat('0',1,nZerosToAdd),num2str(n)];
end