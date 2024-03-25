function [TS2, TS2pos] = match_TS(TS1,spots1,spots2,TS_radius,zxratio)
    %TS1 holds the indices of the rows corresponding to transcription site
    %wtithin the array TS1
    %The program looks within the array TS2 and selects the spots in the
    %vicinity of each txn soite from TS1
    
    %TS2: cell array; each entry contains the indices of the spots that
    %belong to the corresponding TS in TS1
    
    %TS2: cell array of the position and intensity of the TSs
    
    %TS radius in voxels sets the volume around a TS where to look in
    %channel 2.
    
    %zx ratio is the dz/dx ratio of the voxels
    
    %if one spot from channel 2 is in the vicinty of 2 separate TSs from
    %channel 1, I assign it to the closest TS.
    
    
    if numel(TS1) == 0
        % no transcription site
        TS2 = [];
        TS2pos = [];
        all_TS2 = [];
        all_TS2clean = [];
    else
        %find channel 2 spots next to each TS from channel 1
        for i=1:numel(TS1)
            TS1pos{i} = spots1(TS1{i},:);
            is_next_to_TS = zeros(size(spots2,1),1);
            for j = 1:size(TS1{i},1)
                is_next_to_TS =is_next_to_TS + ((sqrt( (spots2(:,1) - TS1pos{i}(j,1)).^2 +  (spots2(:,2) - TS1pos{i}(j,2)).^2 +  zxratio^2*(spots2(:,3) - TS1pos{i}(j,3)).^2)) < TS_radius);
            end
            is_next_to_TS = logical(is_next_to_TS);
            TS2{i} = find(is_next_to_TS);
            TS2pos{i,1} = spots2(is_next_to_TS ,:);
            %generate a column with the TS ID
            if i==1
                TS2idx = ones(numel(TS2{i}),1);
            else
                TS2idx = [TS2idx;i*ones(numel(TS2{i}),1)];
            end
        end

        %remove doublets
        all_TS2 =[cell2mat(TS2pos),TS2idx];

        [~,I,~] = unique(all_TS2(:,1:size(all_TS2,2)-1), 'rows');
        all_TS2clean = all_TS2(I,:);
     
        hasDuplicates = size(all_TS2clean ,1) < size(all_TS2,1);
        
        if hasDuplicates
            %separate TS positions in one array (all_TS2clean) that contains only non
            %repeated values and the other (dupRowValues) all the repeated values
            ixDupRows = zeros(size(all_TS2,1),1);
            for i=1:size(all_TS2clean,1)
                idxtmp = logical((all_TS2(:,1) == all_TS2clean(i,1)).*(all_TS2(:,2) == all_TS2clean(i,2)));
                if sum(idxtmp)>1
                    ixDupRows = logical(ixDupRows + idxtmp);   
                end
            end
            all_TS2clean = all_TS2(~ixDupRows,:);
            dupRowValues = all_TS2(ixDupRows,:);
            
            %loop through the array containing the repeated values
            %for each set of repeated values, assign the spot to the closest TS
            clear_duplicates = [];
            while numel(dupRowValues) > 0
                idxtmp = logical((dupRowValues(:,1) == dupRowValues(1,1)).*(dupRowValues(:,2) == dupRowValues(1,2)))';
                duplicate_set = dupRowValues(idxtmp,:);
                rmin = TS_radius+1;
                imin = 1;
                for i=1:size(duplicate_set,1)
                    rtmp = sqrt( (duplicate_set(i,1) - TS1pos{duplicate_set(i,end)}(1,1)).^2 +...
                        (duplicate_set(i,2) - TS1pos{duplicate_set(i,end)}(1,2)).^2 + ...
                        (zxratio)^2*(duplicate_set(i,3) - TS1pos{duplicate_set(i,end)}(1,3)).^2) < TS_radius;
                    if rtmp<rmin
                        imin = i;
                    end
                end
                if numel(clear_duplicates) == 0
                    clear_duplicates = duplicate_set(imin,:);
                else
                    clear_duplicates = [clear_duplicates;duplicate_set(imin,:)];
                end
                
                dupRowValues = dupRowValues(~idxtmp,:);
                duplicate_set = duplicate_set(duplicate_set(:,end) ~= duplicate_set(:,1));
            end
            all_TS2clean = [all_TS2clean;clear_duplicates];
            all_TS2clean = sortrows(all_TS2clean,8);
            
            %clean up the TS and TSpos cell arrays
            TS2 = []; TS2pos = [];
            for i=1:numel(TS1)
                TS2{i,1} = find(all_TS2clean(:,8) == i);
                TS2pos{i,1} = all_TS2clean(TS2{i,1}(:) ,:);
            end
            
        end
    end

end