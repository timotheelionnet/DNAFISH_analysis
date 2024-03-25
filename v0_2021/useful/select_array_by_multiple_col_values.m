function sel_arr = select_array_by_multiple_col_values(arr,range)
%inclusive
%selects the rows of an array according to threshold value on the values
%of the different columns.
%arr is the array, 
%range is the array formatted as:
%   [ col1_min , col2_min , ... ;
%     col1_max , col2_max , ...]

ny = size(arr,2);
ncol = size(range,2);
if ny<ncol
    disp('your array does not have enough columns - get yourslef a greek temple');
end

sel_arr = arr;
for i =1:ncol
    sel_arr = select_array_by_col_value(sel_arr,i,range(1,i),'min');
    if sel_arr == 0
        return
    end
    sel_arr = select_array_by_col_value(sel_arr,i,range(2,i),'max');  
    if sel_arr == 0
        return
    end
end

end

