function loc = sort_loci_into_cell_masks(loc,refImg)
    %take 2D or 3D loc3 and sort spots into 2D cell masks
    
    %check if loc is empty
    if ~isempty(loc)
        % into ROIs as NaN in a newly created column
        outOfBoundsIdx = ...
            loc(:,1) <= 0 | loc(:,1) > size(refImg,1) | ...
            loc(:,2) <= 0 | loc(:,2) > size(refImg,2);
        loc(outOfBoundsIdx,end+1) = NaN;
                
        linIdx = sub2ind(size(refImg),...
            ceil(loc(~outOfBoundsIdx,1)),...
            ceil(loc(~outOfBoundsIdx,2)));
                
        % assign loci to cells/nuclei
        loc(~outOfBoundsIdx,end) = refImg(linIdx);
    else
        disp('loc3 file is empty!')
    end
end