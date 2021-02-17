function n = find_last_index_above_value(arr,value)
%inclusive
%assumes the array is ordered with decreasing values
nmin = 1;
nmax = size(arr,1);

if max(arr(:,1)) < value
    n = nmax;
    return;
end

while nmax-nmin>1
    nmed = round((nmin+nmax)/2);
    if arr(nmed,1) < value
        nmax = nmed;
    elseif arr(nmed,1) >= value
        nmin = nmed;
    end
end

if arr(nmax) <= value
    n = nmin;
else
    n = nmax;
end

end