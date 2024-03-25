addpath('../fileSet');
addpath('../iniconfig');
addpath('../cbrewer');

%% get list of files
try
    configFileName;
catch
    % config file path
    configFileName = '20210302_DNA_FISH_cell_cycle_config.ini';
end
% create file list object
fs = fileSet;
fs.buildFileSetFromConfig(configFileName);

% reformat channel, replicate and FOV as numbers
convertNonNumeralStrings = 1; 
fs.toNumeric({'Channel','TechnicalReplicate','FOV'},convertNonNumeralStrings);

%read parameters from config file relative to analysis into a parameters object
params = readParamsFromIniFile(configFileName);

%% sort cells to G1/2 phase by DNA FISH spots
curFileList = fs.getFileName({},{},'match','all');

for i = 1:length(curFileList)
    %load match file
    match = load(curFileList{i},'-ascii');
    [~,f,~] = fileparts(curFileList{i});
    if isempty(match)
        continue
    end
    %loop through all cells
    cellCycle = nan(length(match(:,7)),1);
    for j = 1:length(match(:,7))
        count = length(find(match(:,7) == match(j,7)));
        if count <= 2
            cellCycle(j) = 1;
        else
            if count<= 4
                cellCycle(j) = 2;
            end
        end
    end
    match = [match,cellCycle];
    
    if ~ exist(fullfile(params.Output.outFolder,'loci_match_cell_cycle'),'dir')
        mkdir(fullfile(params.Output.outFolder,'loci_match_cell_cycle'));
    end
    matchFileName = fullfile(params.Output.outFolder,'loci_match_cell_cycle',[f,'.dloc3']);
    save(matchFileName,'match','-ascii');
end


%% plot DAPI/EU intensity grouped by cell cycle
I = readmatrix(fullfile(params.InputFolders.MaskFolder,'Cells.csv'));
curFileList = fs.getFileName({},{},'match','all');

% For CUTLL
DAPI1 =[];
DAPI2 =[];
EU1 = [];
EU2 =[];
count1 = 0;
count2 = 0;
count = 0;
for i = 1:length(I(:,1))
    %load match file
    [~,f,~] = fileparts(curFileList{I(i,1)});
    if contains(f,'CUTLL')
        matchFileName = fullfile(params.Output.outFolder,'loci_match_cell_cycle',[f,'.dloc3']);
        try
            match = load(matchFileName,'-ascii');
        catch
            continue
        end
        idx = find(match(:,7)==I(i,2),1);
        
        if match(idx,11) ==1
            DAPI1 = [DAPI1;I(i,13)];
            EU1 = [EU1;I(i,12)];
            count1 =count1+1;
        else
            if match(idx,11) ==2
                DAPI2 = [DAPI2;I(i,13)];
                EU2 = [EU2;I(i,12)];
                count2= count2+1;
            end
        end
        count = count +1;
    end
end

fprintf('There are %d CUTLL cells in total, %d are in G1 phase, %d are in G2 phase.\n',count,count1,count2)

%plot
plot_fig({DAPI1,DAPI2},{'G1','G2'},'Integrated Intensity',[],50,params.Output.outFolder,'DAPI_CUTLL_total')
plot_fig({EU1,EU2},{'G1','G2'},'Integrated Intensity',[],50,params.Output.outFolder,'EU_CUTLL_total')

% For T cells
DAPI1 =[];
DAPI2 =[];
EU1 = [];
EU2 =[];
count1 = 0;
count2 = 0;
count = 0;
for i = 1:length(I(:,1))
    %load match file
    [~,f,~] = fileparts(curFileList{I(i,1)});
    if ~contains(f,'CUTLL')
        matchFileName = fullfile(params.Output.outFolder,'loci_match_cell_cycle',[f,'.dloc3']);
        try
            match = load(matchFileName,'-ascii');
        catch
            continue
        end
        idx = find(match(:,7)==I(i,2),1);

        if match(idx,11) ==1
            DAPI1 = [DAPI1;I(i,13)];
            EU1 = [EU1;I(i,12)];
            count1 =count1+1;
        else
            if match(idx,11) ==2
                DAPI2 = [DAPI2;I(i,13)];
                EU2 = [EU2;I(i,12)];
                count2= count2+1;
            end
        end
        count = count +1;
    end
end

fprintf('There are %d T cells in total, %d are in G1 phase, %d are in G2 phase.\n',count,count1,count2)

%plot
plot_fig({DAPI1,DAPI2},{'G1','G2'},'Integrated Intensity',[],50,params.Output.outFolder,'DAPI_T_cells_total')
plot_fig({EU1,EU2},{'G1','G2'},'Integrated Intensity',[],50,params.Output.outFolder,'EU_T_cells_total')

%% plot CUTLL centroid distances grouped by treatments and cell cycle
try
    cond;
    colorID;
catch
    %find list of conditions
    cond = find_cond(fs.fList);
    colorID = [];
    if params.channelDescription.ConditionDescription
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
end

%read in centroid dis
disc1 = {};
disc2 = {};
for k = 1:numel(cond(:,1))
    curFileList = fs.getFileName({'Condition'},cond(k),'match','all');
    count1 = 0;
    count2 = 0;
    count = 0;
    dis1 = [];
    dis2 = [];
    for j = 1:length(curFileList)
        %load distance
        [d,f,~] = fileparts(curFileList{j});
        try
            dist = load(fullfile(d,strcat(f,'_cen_dist.txt')),'-ascii');
        catch
            continue
        end        
        %load match file
        matchFileName = fullfile(params.Output.outFolder,'loci_match_cell_cycle',[f,'.dloc3']);
        try
            match = load(matchFileName,'-ascii');
        catch
            continue
        end
        for i = 1:length(match(:,1))
            if match(i,11) ==1
                dis1 = [dis1;dist(i)];
                count1= count1+1;
            else
                if match(i,11) ==2
                    dis2 = [dis2;dist(i)];
                    count2= count2+1;
                end
            end
            count = count +1;
        end
    end
    disc1{k}=dis1;
    disc2{k}=dis2;
    fprintf('%s: %d cells in total, %d are in G1 phase, %d are in G2 phase.\n',cond{k},count,count1,count2)
end

plot_fig(disc1,cond,'distance/nm',colorID,params.Settings.hist_binSize,params.Output.outFolder,'dis_cen_G1')
plot_fig(disc2,cond,'distance/nm',colorID,params.Settings.hist_binSize,params.Output.outFolder,'dis_cen_G2')

%% functions for plot
function cond = find_cond(fList)
    %find list of conditions
    cond = unique(fList(:,3));
    cond = cond.(1);%convert table to array
end
function cond2=group_conditions(cond,key)
    cond1 = strings(length(cond),length(key));
    k = zeros(1,length(key));
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
    
    %sort by drug concentration
    cond2 = strings(size(cond1));
    for i = 1:length(key)
        con = regexp(cond1(:,i),'\d+(?:\.\d+)?uM','match');
        con = double(erase([con{:}],'uM'));
        [~,idx] = sort(con);
        for j = 1:k(i)
            if contains(cond1(j,i),'ctrl')
                cond2(1,i)=cond1(j,i);
            else
                cond2(find(idx==j,1)+1,i)=cond1(j,i);
            end
        end
    end
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
    array2 = array2((~strcmp(array,""))&(~ismissing(array)));
    array = array((~strcmp(array,""))&(~ismissing(array)));%remove "" and <missing>
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
            h = plot(x,n);
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
    try
        legend(cond)
        xlabel(label)
    catch
    end
    ylabel('Probability')
    hold off
    savefig(figure(1),fullfile(outFolder,['PDF_',figname,'.fig']))
        
    %save CDF
    figure(2)
    try
        legend(cond,'Location','southeast')
        xlabel(label)
    catch
    end
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
