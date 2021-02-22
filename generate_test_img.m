saving_folder = 'R:\lionnt01lab\lionnt01labspace\Yi_Fu\test'; 

ims = 100; %generate 100 images
imSize3D = [500,500,20];  % size in pixel of image
voxSize3D = [100,100,500]; % size in nm of pixels in x,y,z (can be different)
nPairs3D = 20; % number of spot pairs to fit in the image
meanSeparationDist3D = 0; % average distance between two spots in a pair in nm
stdSeparationDist3D = 100; % std of the distance between two spots in a pair in nm
meanSpotSize3D = [200,500,1000]; % average spot size in nm in x,y,z (can be different)
stdSpotSize3D = [200,200,200]; % std of spot size in nm in x,y,z (can be different)
meanSpotInt3D = 1000; % average spot intensity (amplitude of the 2D gaussian)
stdSpotInt3D = 300; % std of the spot intensity
noiseStd3D = 100; % std of gaussian noise added to image

for i = 1:ims
    [img3D1,img3D2,spotPositions3D,spotSizes3D,spotIntensitie3Ds,spotDensity3D] = ...
        generateSpotPairTestImg(imSize3D,voxSize3D,nPairs3D,...
        meanSeparationDist3D,stdSeparationDist3D,meanSpotInt3D,stdSpotInt3D,...
        meanSpotSize3D,stdSpotSize3D,noiseStd3D);

    % compute actual distances between centroids of generated spots
    x3D = sqrt( (spotPositions3D(:,1) - spotPositions3D(:,4)).^2 + ...
        (spotPositions3D(:,2) - spotPositions3D(:,5)).^2 + ...
        (spotPositions3D(:,3) - spotPositions3D(:,6)).^2);

    % plot 2D projection of spot positions in nm
    figure('name','3D spot pairs - positions in nm');
    hold;
    scatter(spotPositions3D(:,1),spotPositions3D(:,2));
    scatter(spotPositions3D(:,4),spotPositions3D(:,5));

    % save images
    save_as_tiff(img3D1,fullfile(saving_folder,strcat('C1-',string(i),'.tif')));
    save_as_tiff(img3D2,fullfile(saving_folder,strcat('C2-',string(i),'.tif')));

    % save x3D
    save(fullfile(saving_folder,strcat('centroid_',string(i),'.txt')),'x3D','-ascii')
end