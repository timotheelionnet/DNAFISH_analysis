

%%
warning('off','all')
outFolder = fullfile(myParams.Output.outFolder, "res", "sorted_spots");

% loop thru data by FOV

my_FOVs = unique(fs.fList.FOV);


for i = 1:numel(my_FOVs)
% for i = 13:13

    cur_FOV = my_FOVs(i);

    % loop thru each channel for current FOV 
    for j = 1:numel(allFISHChannels)
%     for j =1:1

        cur_ch = allFISHChannels(j);

        % get loc file

        cur_loc_path = string(fs.getFileName({'FOV', 'Channel'}, {cur_FOV, cur_ch}, {file2filter}));
        cur_loc = load(string(fs.getFileName({'FOV', 'Channel'}, {cur_FOV, cur_ch}, {file2filter})));

        % make a copy for downstream filtering
        cur_loc_new = cur_loc;
        
        % get mask
        cur_mask_path = fs.getFileName({'FOV', 'Channel'}, {cur_FOV, cur_ch}, 'mask');
        cur_mask_img = readTifStackWithImRead(cur_mask_path{1});
        
        % get click n keep table
        cur_cnk = readtable(string(fs.getFileName({'FOV', 'Channel'}, {cur_FOV, cur_ch}, {'cNk'})));

        % check if empty
        col_names = string(cur_cnk.Properties.VariableNames);

        if ~contains(col_names, 'X') % empty tables have no X variable
            keep_coords = [];

        
        else
            
            % Y and X are swapped versus FIJI
       
            keep_coords = [ cur_cnk.Y, cur_cnk.X ];

        end
        
       


        % column index number where cellIDs are stored
        cellID_idx = width(cur_loc) + 1;
        
        % sort the loci
        sorted_loc = sort_loci_into_cell_masks(cur_loc, cur_mask_img);

        % firure out cellIDs of the cells we clicked on
        cells2keep = sort_loci_into_cell_masks(keep_coords, cur_mask_img); %% the cellID will be in col 3
        

        if ~isempty(cells2keep) & ~isempty(sorted_loc)
            % remove loc that are not in desired cells
        
            rmv_idx = ~(ismember(sorted_loc(:, cellID_idx), cells2keep(:, 3)));
    
            cur_loc_new(rmv_idx, :) = [];
        else
            % there are no cells to keep
            cur_loc_new(:,:) = [];

        end

        [f,n,x] = fileparts(cur_loc_path{1});
        saveName = fullfile(outFolder, strcat(n, x));
        save(saveName,'cur_loc_new', '-ascii');

        % disp('saved')
        % saveName
        % cur_loc_new

         if ~isempty(cur_loc) & ~isempty(cur_loc_new) & j ==1 & mod(i,10) == 0 % lets look at every 10th FOV
            % visualize
            figure
            
            imagesc(cur_mask_img)
            hold on

            colormap("colorcube")
          
            
            scatter( keep_coords(:,2), keep_coords(:,1), 'hexagram', "MarkerFaceColor", 'white', "MarkerEdgeColor", "black", DisplayName="clicked" )
            hold on 
            scatter(cur_loc(:,2), cur_loc(:,1), 'filled', "diamond", "MarkerFaceColor", "red", "MarkerEdgeColor","black", DisplayName="discard")
            hold on
            scatter( cur_loc_new(:,2), cur_loc_new(:,1), 'diamond', "MarkerFaceColor", 'green', "MarkerEdgeColor", "black", DisplayName="kept" )
            
            
            hold on

            
            legend('Location', 'northeast')
            title('ClickNKeep Cells--FOV ', string(i))
        end


        % find where we are in table and add reg loc3 to fl
        % first need to find col idx of the file type we are sorting
        curFileType_idx = find(ismember(fs.fList.Properties.VariableNames, file2sort));
        % then we find what row we are in
        tempFl = fs.fList(:,curFileType_idx);
        curFl_row_idx = find(ismember(tempFl{:,1}, curLocFl{j}));
        
        addFileFl_idx =  find(ismember(fs.fList.Properties.VariableNames, string(file2sort)));
        % ADD function to SAVE
        fs.fList(curFl_row_idx, addFileFl_idx) = {char(saveName)};

        





       

        

    end

end

warning('on','all')