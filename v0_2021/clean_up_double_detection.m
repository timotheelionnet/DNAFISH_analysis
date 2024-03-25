function cleanSpots = clean_up_double_detection(loc3nm, threshold, nDims)

    [dr, ~] = compute_distMatrix_from_loc(loc3nm,[],nDims);
    toKeep = [];
    for i=1:size(dr,1)
        spotsWithinThreshold = find( dr(i,:) < threshold );
        if isempty(spotsWithinThreshold)
            if ~ismember(i,toKeep)
                toKeep = [toKeep;i];
            end
        else
            intensities = loc3nm(spotsWithinThreshold, nDims+1);
            [~,j] = max(intensities);
            if ~ismember(j,toKeep)
                toKeep = [toKeep;j];
            end
        end
    end

    cleanSpots = loc3nm(toKeep,:);

end