function odata = bin_array(idata,npts,mode)

nbins = ceil(size(idata,1)/npts);


if strcmp(mode,'mean')
    for i=1:nbins-1
        odata(i,1:size(idata,2)) = mean(idata((i-1)*npts+1:i*npts,:));    
        odata(i,size(idata,2)+1:2*size(idata,2)) = std(idata((i-1)*npts+1:i*npts,:));    
    end
    odata(nbins,1:size(idata,2)) = mean(idata((nbins-1)*npts+1:size(idata,1),:));
    odata(nbins,size(idata,2)+1:2*size(idata,2)) = std(idata((nbins-1)*npts+1:size(idata,1),:)); 

elseif strcmp(mode,'median')
    for i=1:nbins-1
        odata(i,1:size(idata,2)) = median(idata((i-1)*npts+1:i*npts,:));    
        odata(i,size(idata,2)+1:2*size(idata,2)) = std(idata((i-1)*npts+1:i*npts,:));    
    end
    odata(nbins,1:size(idata,2)) = median(idata((nbins-1)*npts+1:size(idata,1),:));
    odata(nbins,size(idata,2)+1:2*size(idata,2)) = std(idata((nbins-1)*npts+1:size(idata,1),:)); 
end

end