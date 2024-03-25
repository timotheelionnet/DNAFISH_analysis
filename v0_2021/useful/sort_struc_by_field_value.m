function struc2 = sort_struc_by_field_value(struc1,fieldstring,orientation)

field_arr = eval(['[struc1.',fieldstring,']'])';

[sort_arr,indices] = sort_array_by_col_value(field_arr,1,orientation);

struc2 = struc1(indices);

clear('sort_arr','field_arr','indices');
end