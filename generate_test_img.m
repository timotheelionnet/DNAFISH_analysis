saving_folder = 'Z:\lionnt01lab\lionnt01labspace\Yi_Fu\test\042021SimulateDNAFISH\400'; 

ims = 100; %generate 100 images
imSize3D = [500,500,20];  % size in pixel of image
voxSize3D = [73,73,250]; % size in nm of pixels in x,y,z (can be different)
nPairs3D = 20; % number of spot pairs to fit in the image
meanSeparationDist3D = 400; % average distance between two spots in a pair in nm
stdSeparationDist3D = 700; % std of the distance between two spots in a pair in nm
meanSpotSize3D = [150,150,750]; % average spot size in nm in x,y,z (can be different)
stdSpotSize3D = [30,30,25]; % std of spot size in nm in x,y,z (can be different)
meanSpotInt3D = 12000; % average spot intensity (amplitude of the 2D gaussian)
stdSpotInt3D = 500; % std of the spot intensity
noiseStd3D = 500; % std of gaussian noise added to image

for i = 1:ims
    [img3D1,img3D2,spotPositions3D,spotSizes3D,spotIntensitie3Ds,spotDensity3D,imMask] = ...
        generateSpotPairTestImg2(imSize3D,voxSize3D,nPairs3D,...
        meanSeparationDist3D,stdSeparationDist3D,meanSpotInt3D,stdSpotInt3D,...
        meanSpotSize3D,stdSpotSize3D,noiseStd3D);

    % compute actual distances in nm between centroids of generated spots
    x3D = sqrt( voxSize3D(1)^2*(spotPositions3D(:,1) - spotPositions3D(:,4)).^2 + ...
        voxSize3D(2)^2*(spotPositions3D(:,2) - spotPositions3D(:,5)).^2 + ...
        voxSize3D(3)^2*(spotPositions3D(:,3) - spotPositions3D(:,6)).^2);

    % save images
    save_as_tiff(imMask,fullfile(saving_folder,strcat('C1-',string(i),'max_Mask.tiff')));
    save_as_tiff(img3D1,fullfile(saving_folder,strcat('C2-',string(i),'.tif')));
    save_as_tiff(img3D2,fullfile(saving_folder,strcat('C3-',string(i),'.tif')));
    
    pairs = [spotPositions3D,x3D];

    % save x3D
    save(fullfile(saving_folder,strcat('centroid_',string(i),'.txt')),'pairs','-ascii')
end