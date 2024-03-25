function sel_arr = select_array_by_col_value(arr,ncol,value,orientation)

%selects the rows of an array according to a threshold value on the value
%sof 1 column.
%arr is the array, 
%col is the column which values are selected from
%orientation is a string, either 'min', 'max' or 'equal'. It indicates
%whether the value is the resp. the min/max/exact match of the selected row
%values
%comparison is inclusive

ny = size(arr,2);
if ny<ncol
    disp('your array does not have enough columns - get yourslef a greek temple');
end

if strcmp(orientation,'min')
    [sortVal indices] = sort(arr(:,ncol),'descend');
    nx = size(arr,1);
    sort_arr(1:nx,1:ny) = arr(indices(1:nx),1:ny);
    if sort_arr(1,ncol)<value
        sel_arr = 0;
        return
    elseif sort_arr(nx,ncol)>value
        sel_arr = sort_arr;
        return
    end
    n = find_last_index_above_value(sort_arr(:,ncol),value); %dichotomy search
    sel_arr = sort_arr(1:n,1:ny);
    
elseif strcmp(orientation,'max')
    [sortVal indices] = sort(arr(:,ncol),'ascend');
    nx = size(arr,1);
    sort_arr(1:nx,1:ny) = arr(indices(1:nx),1:ny);
    if sort_arr(1,ncol)>value
        sel_arr = 0;
        return
    elseif sort_arr(nx,ncol)<value
        sel_arr = sort_arr;
        return
    end
    n = find_last_index_below_value(sort_arr(:,ncol),value); %dichotomy search
    sel_arr = sort_arr(1:n,1:ny);
    
elseif strcmp(orientation,'equal')
    nx = size(arr,1);
    sel_arr = zeros(1,ny);
    j=1;
    for i= 1:nx
        if arr(i,ncol)==value
            sel_arr(j,:) = arr(i,:);
            j = j+1;
        end
    end
    
end

clear ('sortVal','sort_arr');
end

