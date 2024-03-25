function n = find_last_index_below_value_exclusive(arr,value)
%assumes the array is ordered with increasing values
%exclusive (finds last index < value)
%returns 0 if all array entries are greater than value

nmin = 1;
nmax = size(arr,1);

if min(arr(:,1)) > value
    n= 0;
    return;
end

while nmax - nmin>1
    nmed = round((nmin+nmax)/2);
    if arr(nmed,1) >= value
        nmax = nmed;
    elseif arr(nmed,1) < value
        nmin = nmed;
    end
end

if arr(nmax) >= value
    n = nmin;
else
    n = nmax;
end

end