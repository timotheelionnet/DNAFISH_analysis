function crop_box(ims, loc, outDir, boxSize, thickness)
    %crop a box from ims using loc location, discard boxes that is not
    %fully in the image
    %Output cropped boxes in folders
    %Output a new loc file contain only loci fully in the ims
    
    locin = [];
    loci = cell('');
    j = 0;
    nDigitsMax = ceil(log(size(loc,1)+1)/log(10));    
    wx = boxSize(1);
    wy = boxSize(2);
    wz = boxSize(3);
    
    for i=1:size(loc,1)
    
        is_fully_in_image = (loc(i,1) >= 1+wx) && ...
            (loc(i,1) <= size(ims,1) - wx) && ...
            (loc(i,2) >= 1+wy) && ...
            (loc(i,2) <= size(ims,2) - wy) && ...
            (loc(i,3) >= 1+wz) && ...
            (loc(i,3) <= size(ims,3) - wz);
        
        if is_fully_in_image
            j = j+1;
            locin = [locin; loc(i,:)];
            loci{j} = collect_box(ims,loc(i,1:3),[wx,wy,wz]);
            nDigits = ceil(log(j+1)/log(10));
            nZerosToAdd = nDigitsMax - nDigits;
            numString = [repmat('0',1,nZerosToAdd),num2str(j)];

            [~,stack_bg_corr] = ...
                gen_linear_interpol_clean_smallYi(ims,loc(i,1:3),[wx,wz],thickness,'small');
            save_as_tiff(loci{j},fullfile(outDir,['locus',numString,'.tif']));
            save_as_tiff(stack_bg_corr,fullfile(outDir,['locus_bgcorr',numString,'.tif']));
        end
        
    end

    locFileName = fullfile(outDir,'location.loc3');
    save(locFileName,'locin','-ascii');
    disp(['collected ',num2str(j),' entire bead substacks from image; '...
        'eliminated ',num2str(size(loc,1)-j),' beads near stack edges.']);
end