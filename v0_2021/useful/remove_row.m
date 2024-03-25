function out = remove_row(arr,row_num)
    nr = size(arr,1);
    arr(row_num:nr-1,:) = arr(row_num+1:nr,:);
    out = arr(1:nr-1,:);
end