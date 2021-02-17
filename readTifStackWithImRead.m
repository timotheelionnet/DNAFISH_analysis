function tiff_stack = readTifStackWithImRead(fileName)
tiff_info = imfinfo(fileName); % return tiff structure, one element per image
tiff_stack = zeros(tiff_info(1).Width,tiff_info(1).Height,size(tiff_info, 1)); % read in first image
%concatenate each successive tiff to tiff_stack
for i = 1 : size(tiff_info, 1)
    tiff_stack(:,:,i) = imread(fileName, i);
end



