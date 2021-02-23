addpath('../fileSet');
addpath('../iniconfig');
addpath('../cbrewer');

%% get list of files

% config file path
configFileName = '20210203_DNA_FISH_plot_config.ini';

% create file list object
fs = fileSet;
fs.buildFileSetFromConfig(configFileName);

% reformat channel, replicate and FOV as numbers
convertNonNumeralStrings = 1; 
fs.toNumeric({'Channel','TechnicalReplicate','FOV'},convertNonNumeralStrings);

%read parameters from config file relative to analysis into a parameters object
params = readParamsFromIniFile(configFileName);

%% Plot PDF and CDF for dis between airlocalize center
%find list of conditions
cond = find_cond(fs.fList);

%group conditions and generate color scheme for ploting
des = parse_description(params.channelDescription.ConditionDescription);
cond = group_conditions(cond,des);
size0 = size(cond);
color1 = color_generator(size0);
idx = 1:numel(cond);
idx = reshape(idx,size0);
[cond,idx1] = flat_array(cond,'column',string(idx));
[row,col]=ind2sub(size0,str2double(idx1));
colorID = zeros(length(idx1),3);
for i = 1:length(idx1)
    colorID(i,:) = color1(row(i),col(i),:);
end
cond = mat2cell(cond,ones(length(cond),1));

%read in dis
dist = {};
for i = 1:numel(cond(:,1))
    curFileList = fs.getFileName({'Condition'},cond(i),'match','all');
    dis = [];
    for j = 1:length(curFileList)
        %load match file
        match = load(curFileList{j},'-ascii');
    
        %save dist
        dis = [dis;match(:,6)];    
    end
    dist{i} = dis;%each column is dis of a different condition
end

plot_fig(dist,cond,'distance/nm',colorID,params.Settings.hist_binSize,params.Output.outFolder,'dis_Airlocalize')

%% Plot PDF and CDF for dx, dy, dz of airlocalize center
voxsize = [params.Settings.voxSize_dxy,...
            params.Settings.voxSize_dxy,...
            params.Settings.voxSize_dz];
        
try
    colorID;
catch
    %find list of conditions
    cond = find_cond(fs.fList);
    des = parse_description(params.channelDescription.ConditionDescription);
    cond = group_conditions(cond,des);
    size0 = size(cond);
    color1 = color_generator(size0);
    idx = 1:numel(cond);
    idx = reshape(idx,size0);
    [cond,idx1] = flat_array(cond,'column',string(idx));
    [row,col]=ind2sub(size0,str2double(idx1));
    colorID = zeros(length(idx1),3);
    for i = 1:length(idx1)
        colorID(i,:) = color1(row(i),col(i),:);
    end
    cond = mat2cell(cond,ones(length(cond),1));
end

%read in dis
dx = {};
dy = {};
dz = {};
for i = 1:numel(cond(:,1))
    curFileList = fs.getFileName({'Condition'},cond(i),'match','all');
    dis = [];
    for j = 1:length(curFileList)
        %load match file
        match = load(curFileList{j},'-ascii');
    
        %save dist
        dis = [dis;match(:,3:5)];    
    end
    dis = int64(convert_loc_pix_to_nm(dis,1./voxsize));
    dx{i} = double(dis(:,1));%each column is dis of a different condition
    dy{i} = double(dis(:,2));
    dz{i} = double(dis(:,3));
end

plot_fig(dx,cond,'pixel',colorID,1,params.Output.outFolder,'dx')
plot_fig(dy,cond,'pixel',colorID,1,params.Output.outFolder,'dy')
plot_fig(dz,cond,'pixel',colorID,1,params.Output.outFolder,'dz')

%% plot PDF and CDF for overlap

try
    colorID;
catch
    %find list of conditions
    cond = find_cond(fs.fList);
    des = parse_description(params.channelDescription.ConditionDescription);
    cond = group_conditions(cond,des);
    size0 = size(cond);
    color1 = color_generator(size0);
    idx = 1:numel(cond);
    idx = reshape(idx,size0);
    [cond,idx1] = flat_array(cond,'column',string(idx));
    [row,col]=ind2sub(size0,str2double(idx1));
    colorID = zeros(length(idx1),3);
    for i = 1:length(idx1)
        colorID(i,:) = color1(row(i),col(i),:);
    end
    cond = mat2cell(cond,ones(length(cond),1));
end

%read in overlap
dist = {};
for i = 1:numel(cond(:,1))
    curFileList = fs.getFileName({'Condition'},cond(i),'match','all');
    dis = [];
    for j = 1:length(curFileList)
        %load overlap file
        [~,f,~] = fileparts(curFileList{j});
        overlap = load(fullfile(params.InputFolders.OverlapFolder,[f,'_overlap.txt']),'-ascii');
    
        %save dist
        dis = [dis;round(overlap,10)];    
    end
    dis0 = 1./dis;
    dis0(isinf(dis0)) = 100000;
    dis0(dis0>100000) = 100000;
    dist{i} = dis0;%each column is dis of a different condition
end

plot_fig(dist,cond,'1/correlation',colorID,100,params.Output.outFolder,'corr')

%% Plot PDF and CDF for dis between loci center
voxsize = [params.Settings.voxSize_dxy,...
            params.Settings.voxSize_dxy,...
            params.Settings.voxSize_dz];

try
    colorID;
catch
    %find list of conditions
    cond = find_cond(fs.fList);
    des = parse_description(params.channelDescription.ConditionDescription);
    cond = group_conditions(cond,des);
    size0 = size(cond);
    color1 = color_generator(size0);
    idx = 1:numel(cond);
    idx = reshape(idx,size0);
    [cond,idx1] = flat_array(cond,'column',string(idx));
    [row,col]=ind2sub(size0,str2double(idx1));
    colorID = zeros(length(idx1),3);
    for i = 1:length(idx1)
        colorID(i,:) = color1(row(i),col(i),:);
    end
    cond = mat2cell(cond,ones(length(cond),1));
end

%calculate dis between center
dist = {};
distx = {};
disty = {};
distz = {};
for i = 1:numel(cond(:,1))
    curFileList = fs.getFileName({'Condition'},cond(i),'match','all');
    dis1 = [];
    dis2 = [];
    for j = 1:length(curFileList)
        %load match file
        match = load(curFileList{j},'-ascii');        
        [d,f,~] = fileparts(curFileList{j});
        disp('processing file:')
        disp(f)
        
        % find the conditions in the current file
        fl = fs.getFileName({},{},'match','all');
        idx = find( ismember(fl,curFileList{j}) );
        curCond = {};
        for k = 1:numel(fs.conditions)
            curCond = [curCond,fs.fList.(fs.conditions{k})(idx)];
        end
        %load loc3 of centroid
        C1 = curCond{1};
        C2 = curCond{2};
        CropFolder1 = fullfile(params.InputFolders.CropFolder,...
            strcat('C',C1,f(5:length(f)),'_crop'));
        CropFolder2 = fullfile(params.InputFolders.CropFolder,...
            strcat('C',C2,f(5:length(f)),'_crop'));
        cen1 = load(fullfile(CropFolder1,'centroid.loc3'),'-ascii');
        cen2 = load(fullfile(CropFolder2,'centroid.loc3'),'-ascii');
        deltanm = match(:,3:5);
        delta = int64(convert_loc_pix_to_nm(deltanm,1./voxsize));
        
        dis = [];
        
        for k = 1:length(match)
            dis3 = cen2(match(k,2),:)-cen1(match(k,1),:)+double(delta(k,:));
            dis2 = [dis2;dis3];
            dis3 = convert_loc_pix_to_nm(dis3,voxsize);
            dis0 = sqrt(sum(dis3.^2));
            %save dist
            dis = [dis;dis0];
        end
        save(fullfile(d,strcat(f,'_cen_dist.txt')),'dis','-ascii')
        dis1 = [dis1;dis];
    end
    dist{i} = dis1;%each column is dis of a different condition
    distx{i} = dis2(:,1);
    disty{i} = dis2(:,2);
    distz{i} = dis2(:,3);
end

%plot CDF and PDF
plot_fig(dist,cond,'distance/nm',colorID,params.Settings.hist_binSize,params.Output.outFolder,'dis_cen')
plot_fig(distx,cond,'pixel',colorID,1,params.Output.outFolder,'dx_cen')
plot_fig(disty,cond,'pixel',colorID,1,params.Output.outFolder,'dy_cen')
plot_fig(distz,cond,'pixel',colorID,1,params.Output.outFolder,'dz_cen')

%% Plot PDF and CDF for Rg and volumn
try
    colorID;
catch
    %find list of conditions
    cond = find_cond(fs.fList);
    des = parse_description(params.channelDescription.ConditionDescription);
    cond = group_conditions(cond,des);
    size0 = size(cond);
    color1 = color_generator(size0);
    idx = 1:numel(cond);
    idx = reshape(idx,size0);
    [cond,idx1] = flat_array(cond,'column',string(idx));
    [row,col]=ind2sub(size0,str2double(idx1));
    colorID = zeros(length(idx1),3);
    for i = 1:length(idx1)
        colorID(i,:) = color1(row(i),col(i),:);
    end
    cond = mat2cell(cond,ones(length(cond),1));
end

R1x = {};
R1y = {};
R1z = {};
R2x = {};
R2y = {};
R2z = {};
V1 = {};
V2 = {};
for i = 1:numel(cond(:,1))
    curFileList = fs.getFileName({'Condition'},cond(i),'match','all');
    r1 = [];
    r2 = [];
    v1 = [];
    v2 = [];
    for j = 1:length(curFileList)
        [~,f,~] = fileparts(curFileList{j});
        % find the conditions in the current file
        fl = fs.getFileName({},{},'match','all');
        idx = find( ismember(fl,curFileList{j}) );
        curCond = {};
        for k = 1:numel(fs.conditions)
            curCond = [curCond,fs.fList.(fs.conditions{k})(idx)];
        end
        %load Rg and volume
        C1 = curCond{1};
        C2 = curCond{2};
        CropFolder1 = fullfile(params.InputFolders.CropFolder,...
            strcat('C',C1,f(5:length(f)),'_crop'));
        CropFolder2 = fullfile(params.InputFolders.CropFolder,...
            strcat('C',C2,f(5:length(f)),'_crop'));
        r1 = [r1;load(fullfile(CropFolder1,'Rg.loc3'),'-ascii')];
        r2 = [r2;load(fullfile(CropFolder2,'Rg.loc3'),'-ascii')];
        v1 = [v1;load(fullfile(CropFolder1,'volumn.txt'),'-ascii')];
        v2 = [v2;load(fullfile(CropFolder2,'volumn.txt'),'-ascii')];
    end
    R1x{i} = r1(:,1);
    R1y{i} = r1(:,2);
    R1z{i} = r1(:,3);
    R2x{i} = r2(:,1);
    R2y{i} = r2(:,2);
    R2z{i} = r2(:,3);
    V1{i} = v1;
    V2{i} = v2;
end

plot_fig(R1x,cond,'Rg/pixel',colorID,0.02,params.Output.outFolder,'Rg_x_127')
plot_fig(R2x,cond,'Rg/pixel',colorID,0.02,params.Output.outFolder,'Rg_x_130')
plot_fig(R1y,cond,'Rg/pixel',colorID,0.02,params.Output.outFolder,'Rg_y_127')
plot_fig(R2y,cond,'Rg/pixel',colorID,0.02,params.Output.outFolder,'Rg_y_130')
plot_fig(R1z,cond,'Rg/pixel',colorID,0.02,params.Output.outFolder,'Rg_z_127')
plot_fig(R2z,cond,'Rg/pixel',colorID,0.02,params.Output.outFolder,'Rg_z_130')
plot_fig(V1,cond,'Volume/pixel',colorID,params.Settings.hist_binSize,params.Output.outFolder,'Volume_127')
plot_fig(V2,cond,'Volume/pixel',colorID,params.Settings.hist_binSize,params.Output.outFolder,'Volume_130')

%% functions for plot
function cond = find_cond(fList)
    %find list of conditions
    cond = unique(fList(:,3));
    cond = cond.(1);%convert table to array
end
function cond1=group_conditions(cond,key)
    cond1 = strings(length(cond),length(key));
    k = zeros(length(key));
    for i = 1:length(cond)
        for j = 1:length(key)
            if contains(cond(i),key(j))
                cond1(k(j)+1,j) = string(cond(i));
                k(j) = k(j) + 1;
            end
        end
    end

    %remove empty strings
    cond1(all(strcmp(cond1,""),2),:) = [];
end
function colorID = color_generator(size)
    %generate an array of colors, each column is a different series of
    %colors(max 6), each row has a different shade from darker to brighter
    if size(2) > 6
        disp('Only support 6 color series!')
        colorID = [];
        return
    end
    colorID = zeros(size(1),size(2),3);
    scheme = ["Blues","Oranges","Greens","Purples","Greys","Reds"];
    for i = 1:size(2)
        colors = cbrewer('seq',scheme(i),size(1)+2);
        %discard 2 most light colors and invert, so that darker colors come first
        colorID(:,i,:) = flip(colors(3:length(colors),:),1);
    end
end
function [array, array2] = flat_array(array,direction,array2)
    %flat string array aong specified direction and remove ""
    %Also manipulate array2 as did array
    %Array2 must be same size as array
    
    try
        if size(array2) ~= size(array)
            disp('Sizes of arrays don''t match!')
            return
        end
    catch
        disp('Dimensions of arrays don''t match!')
        return
    end
    if strcmp(direction,'column')
        array = reshape(array,[],1);
        array2 = reshape(array2,[],1);
     
    else
        if strcmp(direction,'row')
            array = reshape(array,1,[]);
            array2 = reshape(array2,1,[]);
        else
            disp('Direction must be by ''column'' or by ''row''.')
            return
        end        
    end
    array2 = array2(~strcmp(array,""));
    array = array(~strcmp(array,""));%remove "" 
end
function plot_fig(varible,cond,label,colorID,binsize,outFolder,figname)
    %This function group varible by conditions and plot a PDF and a CDF for
    %each condition in two plots (one PDF, one CDF)
    %color must be a matrix of RGB value with each row as a color
    %varible is an array with each column content results from a different
    %condition
    %figname is name of figs WITHOUT extention
    
    try
        close(figure(1))
        close(figure(2))
    catch
    end
    
    for i = 1:numel(cond(:,1))
        %plot PDF
        if sum(cell2mat(varible(:,1)) < 0)
            x = cell2mat(varible(:,i));
            binmax = ceil(max(max(x),-min(x))/binsize)*binsize;
            [n,x] = hist(x,-binmax:binsize:binmax);
        else
            [n,x] = hist(cell2mat(varible(:,i)),0:binsize:max(cell2mat(varible(:,i))));
        end
        
        n = n/length(cell2mat(varible(:,i)));
        figure(1);
        try 
            h = plot(x,n,'color',colorID(i,:));            
        catch
            plot(x,n)
        end
        set(h,'LineWidth',1.0)
        hold on
        
        %plot CDF
        figure(2);
        h = cdfplot(cell2mat(varible(:,i)));
        try 
            set(h,'color',colorID(i,:))             
        catch
        end
        set(h,'LineWidth',1.0)
        hold on
    end
 
    %save PDF
    figure(1)
    legend(cond)
    xlabel(label)
    ylabel('Probability')
    hold off
    savefig(figure(1),fullfile(outFolder,['PDF_',figname,'.fig']))
        
    %save CDF
    figure(2)
    legend(cond,'Location','southeast')
    xlabel(label)
    ylabel('Cumulative Probability')
    hold off
    savefig(figure(2),fullfile(outFolder,['CDF_',figname,'.fig']))
end

%%
function des1 = parse_description(des)
    idx = strfind(des,', ');
    des1 = string([]);
    for i = 1:length(idx)
        if i == 1
            des1 = [des1;des(1:idx(1)-1)];
        else
            des1 = [des1;des(idx(i-1)+2:idx(i)-1)];
            if i == length(idx)
                des1 = [des1;des(idx(i)+2:length(des))];
            end            
        end
    end
end