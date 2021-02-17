function arr = replace_NaN_by_zero(arr)

tmp = arr;
arr = zeros(size(arr));
arr(~isnan(tmp)) = tmp(~isnan(tmp));

end



