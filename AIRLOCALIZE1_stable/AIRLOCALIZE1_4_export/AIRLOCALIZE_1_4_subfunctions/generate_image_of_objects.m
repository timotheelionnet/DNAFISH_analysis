function out1 = generate_image_of_objects(cc,p)
    
    %initializing the image as zero
    if p.numdim == 3
        out1 = zeros(p.nx,p.ny,p.nz);
    elseif p.numdim ==2
        out1 = zeros(p.nx,p.ny);
    end

    %filling in objects, assigning a value that is equal to the number of
    %mRNAs:
    %I/Imed is within 0 and 1.5: 1 mRNA
    %I/Imed is within 1.5 and 2.5: 2 mRNAs
    %I/Imed is within 2.5 and 3.5: 3 mRNAs
    %etc
    
    %rejecting the spurious spots
    singles = (cc.INorm > 0) .* (cc.INorm < 1.50).* (cc.FiltObjIdx ~= 0);
    
    nRNA = max([round(cc.INorm).*(cc.FiltObjIdx ~= 0);singles]);
    
    for i = 1:length(cc.PixelIdxList)
        if nRNA(i) ~= 0
            out1(cc.PixelIdxList{1,i}) = nRNA(i);
        end
    end
    


end