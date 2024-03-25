% example script showing a test case of generateSpotPairTestImg

% enter in the line below a folder name to save the images generated
saving_folder = '/Users/lionnt01/Documents/junk/'; 

%% use this section to test the function in 2D
imSize2D = [500,500]; % size in pixel of image
voxSize2D = [100,100]; % size in nm of pixels in x and y (can be different)
nPairs2D = 20; % number of spot pairs to fit in the image
meanSeparationDist2D = 1000; % average distance between two spots in a pair in nm
stdSeparationDist2D = 100; % std of the distance between two spots in a pair in nm
meanSpotSize2D = [200,500]; % average spot size in nm in x and y (can be different)
stdSpotSize2D = [20,20]; % std of spot size in nm in x and y (can be different)
meanSpotInt2D = 1000; % average spot intensity (amplitude of the 2D gaussian)
stdSpotInt2D = 300; % std of the spot intensity
noiseStd2D = 100; % std of gaussian noise added to image


[img2D1,img2D2,spotPositions2D,spotSizes2D,spotIntensities2D,spotDensity2D] = ...
    generateSpotPairTestImg(imSize2D,voxSize2D,nPairs2D,...
    meanSeparationDist2D,stdSeparationDist2D,meanSpotInt2D,stdSpotInt2D,...
    meanSpotSize2D,stdSpotSize2D,noiseStd2D);

%compute actual distances between centroids of generated spots
x2D = sqrt( (spotPositions2D(:,1) - spotPositions2D(:,3)).^2 + ...
    (spotPositions2D(:,2) - spotPositions2D(:,4)).^2);

% plot spot positions in nm
figure('name','2D spot pairs - positions in nm');
hold;
scatter(spotPositions2D(:,1),spotPositions2D(:,2));
scatter(spotPositions2D(:,3),spotPositions2D(:,4));

% save images
save_as_tiff(img2D1,fullfile(saving_folder,'im2D1.tif'));
save_as_tiff(img2D2,fullfile(saving_folder,'im2D2.tif'));

%% use this section to test the function in  3D
imSize3D = [500,500,20];  % size in pixel of image
voxSize3D = [100,100,500]; % size in nm of pixels in x,y,z (can be different)
nPairs3D = 20; % number of spot pairs to fit in the image
meanSeparationDist3D = 1000; % average distance between two spots in a pair in nm
stdSeparationDist3D = 100; % std of the distance between two spots in a pair in nm
meanSpotSize3D = [200,500,1000]; % average spot size in nm in x,y,z (can be different)
stdSpotSize3D = [200,200,200]; % std of spot size in nm in x,y,z (can be different)
meanSpotInt3D = 1000; % average spot intensity (amplitude of the 2D gaussian)
stdSpotInt3D = 300; % std of the spot intensity
noiseStd3D = 100; % std of gaussian noise added to image

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
save_as_tiff(img3D1,fullfile(saving_folder,'im3D1.tif'));
save_as_tiff(img3D2,fullfile(saving_folder,'im3D2.tif'));