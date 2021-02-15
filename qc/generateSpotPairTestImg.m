function [img1,img2,spotPositions,spotSizes,spotIntensities,spotDensity] = ...
    generateSpotPairTestImg(imSize,voxSize,nPairs,...
    meanSeparationDist,stdSeparationDist,meanSpotInt,stdSpotInt,...
    meanSpotSize,stdSpotSize,noiseStd)
% generate an image with nPairs pairs of spots in it, each spot has a
% gaussian distributed intensity profile.

% INPUT
% imSize: size [x,y] or [x,y,z] of image in voxels 
% voxSize: [dx,dy] or [dx,dy,dz] is the vox size in nanometers. If imSize
    % has two elements, dz will be ignored (if present).
% nPairs: number of spot pairs in image
% meanSeparationDist: mean distance in nm of between the 2 spots in each pair
% stdSeparationDist: std of the distance in nm (distances are gaussian
% distributed; set to zero to make the distance deterministic).
% meanSpotSize: [mean_sx,mean_sy,mean_sz]: mean sigma in nm in each dimension of the spot intensity
% profile.
% stdSpotSize: [std_sx,std_sy,std_sz]: std of sigma in nm in each dimension of the spot intensity
% profile (sport sizes are gaussian
% distributed; set to zero to make the size deterministic).
% noiseStd: standard deviation of the gaussian noise added to the image (set to zero to ignore)

% OUTPUT
% img1 is the image with the first half spot pairs
% img2 is the image with the other half spot pairs
% spotPositions lists the positions of the spot pairs in nm units:
    % [x1 y1 x2 y2] in 2D
    % [x1 y1 z1 x2 y2 z2] in 3D
% spotSizes lists the sizes of all the spots:
    % [sx1 sy1 sx2 sy2] in 2D (nm)
    % [sx1 sy1 sz1 sx2 sy2 sz2] in 3D (nm)
% spotIntensities lists the intensity of all the spots:
    % [I1 I2] Intensity value equals the amplitude of the gaussian
% spotDensity is the approximated density of spots in the image
    % it is unitless (basically volue occupied by all spot pairs divided by
    % the total image volume) and is a bit overestimated. Useful to avoid 
    % crowded conditions.
    
% intensity, distances and spot sizes are all gaussian-distributed around 
% their respective means. An offset of 6 noiseSts is added to the images.

% checking that inputs have compatible dimensionality, set nDims
if numel(imSize) == 2
    nDims = 2;
elseif numel(imSize) == 3
    nDims = 3;
    if numel(voxSize) <3 || numel(meanSpotSize) <3
        disp(['Could not create spot pairs, image is 3D but voxel or spot',...
            'inputs is 2D']);
        [img1,img2,spotPositions,spotSizes,spotIntensities,spotDensity] = ...
            deal([],[],[],[],[],[]);
        return
    end
else
    disp(['Could not create spot pairs, image size should be 2 or 3D, '...
        'instead is ',num2str(numel(imSize)),' D']);
    [img1,img2,spotPositions,spotSizes,spotIntensities,spotDensity] = ...
            deal([],[],[],[],[],[]);
    return
end

% crop z dimension if needed, and convert inputs to row vector format
imSize = reshape(imSize(1:nDims),1,[]); 
voxSize = reshape(voxSize(1:nDims),1,[]); 
meanSpotSize = reshape(meanSpotSize(1:nDims),1,[]); 
    
% compute x,y,z radius of a spot pair (in pixel units)
pairWidth = (meanSeparationDist / 2 + 3 * meanSpotSize)./voxSize;

% generate spot pairs in a central region so they do not get cropped
% pairCenters in pixel units
pairCenters = pairWidth + 1 +...
    rand(nPairs,nDims) .* repmat(imSize - 2 * pairWidth -2, nPairs, 1);

% estimate of spot pair density = "volume of a pair"/"acccessible volume"
% volume of a pair is an overestimate so density is also an overestimate.
% if the density value is above 1, some spot pairs will overlap 
% (there is also some chance of overlap at lower density values)
spotDensity = nPairs * ...
    prod(2*pairWidth,'all')/prod(imSize - 2* pairWidth,'all');

% compute coordinates of individual spots in each pair
% gaussian distributed distances (making sure all are >0)
separationDistances = max(0, meanSeparationDist + ...
    randn(nPairs,1)*stdSeparationDist); 

pairAngles = 2*pi*rand(nPairs,nDims-1);

if nDims == 2
    spotPositions = pairCenters - separationDistances/2 .* ...
        [cos(pairAngles) , sin(pairAngles) ]./repmat(voxSize,nPairs,1);
    spotPositions = [spotPositions , pairCenters + meanSeparationDist/2 * ...
        [cos(pairAngles) , sin(pairAngles) ]./repmat(voxSize,nPairs,1)];
else
    spotPositions = pairCenters - separationDistances/2 .*  ...
        [cos(pairAngles(:,1)) .* sin(pairAngles(:,2)) , ...
        sin(pairAngles(:,1)) .* sin(pairAngles(:,2)) ,...
        cos(pairAngles(:,2))]./repmat(voxSize,nPairs,1);
    spotPositions = [spotPositions , pairCenters + separationDistances/2 .*  ...
        [cos(pairAngles(:,1)) .* sin(pairAngles(:,2)) , ...
        sin(pairAngles(:,1)) .* sin(pairAngles(:,2)) ,...
        cos(pairAngles(:,2))]./repmat(voxSize,nPairs,1)];
end


% gaussian distributed intensities (making sure they are > 0)
spotIntensities = abs( meanSpotInt + ...
    randn(nPairs,2)*stdSpotInt ); 

% gaussian distributed spot Sizes in nm (making sure they are > 0)
spotSizes = abs ( repmat(meanSpotSize,nPairs,2) + ...
    randn(nPairs,2*nDims).*repmat(stdSpotSize,nPairs,2) );

% initialize image with gaussian noise (adding 6 sigmas offset to avoid
% nnegative values)
img1 = noiseStd*randn(imSize) + 6*noiseStd;
img2 = noiseStd*randn(imSize) + 6*noiseStd;

for i=1:nPairs
    
    % bounding box coordinates spot 1
    [bmin,bmax] = boundingBoxCoordinates(imSize,spotPositions(i,1:nDims),...
        spotSizes(i,1:nDims),voxSize,nDims);
    
    % calculate gaussian intensity across box spot 1 & add to img 1
    img1 = addGaussianIntensityToImg(img1,bmin,bmax,spotPositions(i,1:nDims),...
        spotIntensities(i,1),spotSizes(i,1:nDims),voxSize,nDims);
    
    % bounding box coordinates spot 2
    [bmin,bmax] = boundingBoxCoordinates(imSize,spotPositions(i,nDims+1:2*nDims),...
        spotSizes(i,nDims+1:2*nDims),voxSize,nDims);
    
    % calculate gaussian intensity across box spot 2 & add to img 2
    img2 = addGaussianIntensityToImg(img2,bmin,bmax,spotPositions(i,nDims+1:2*nDims),...
        spotIntensities(i,2),spotSizes(i,nDims+1:2*nDims),voxSize,nDims);
    
end

% convert spot positions from pixel to nm
spotPositions = spotPositions .* repmat(voxSize,nPairs,2);

end


function img = addGaussianIntensityToImg(img,bmin,bmax,spotPosition,...
    spotIntensity,spotSize,voxSize,nDims)
    
    if nDims == 2
        % create meshgrid swapping x/y to match image coordinates
        [yy,xx] = ...
            meshgrid(bmin(2):bmax(2),bmin(1):bmax(1));
        
        % compute gaussian intensity
        g = exp(-(xx - spotPosition(1)).^2*voxSize(1)^2/(2*spotSize(1)^2));
        g = g .* ...
            exp(-(yy - spotPosition(2)).^2*voxSize(2)^2/(2*spotSize(2)^2));
        g = g * spotIntensity;
        
        img(bmin(1):bmax(1),bmin(2):bmax(2)) = ...
            img(bmin(1):bmax(1),bmin(2):bmax(2)) + g;
    else
        % create meshgrid swapping x/y to match image coordinates
        [yy,xx,zz] = ...
            meshgrid(bmin(2):bmax(2),bmin(1):bmax(1),bmin(3):bmax(3));
        
        % compute gaussian intensity
        g = exp(-(xx - spotPosition(1)).^2*voxSize(1)^2/(2*spotSize(1)^2));
        g = g .* ...
            exp(-(yy - spotPosition(2)).^2*voxSize(2)^2/(2*spotSize(2)^2));
        g = g .* ...
            exp(-(zz - spotPosition(3)).^2*voxSize(3)^2/(2*spotSize(3)^2));
        g = g * spotIntensity;
        
        img(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3)) = ...
            img(bmin(1):bmax(1),bmin(2):bmax(2),bmin(3):bmax(3)) + g;
    end
    
end

function [bmin,bmax] = boundingBoxCoordinates(...
    imSize,spotPosition,spotSize,voxSize,nDims)
% coordinates of bounding box for spots in pixel units
spotSize = 3* spotSize; % taking 3 times the sigma of the gaussian

    bmin(1) = ceil(spotPosition(1) - spotSize(1)/voxSize(1));
    bmin(1) = max( min(bmin(1),imSize(1)), 1);
    bmin(2) = ceil(spotPosition(2) - spotSize(2)/voxSize(2));
    bmin(2) = max( min(bmin(2),imSize(2)), 1);
    bmax(1) = ceil(spotPosition(1) + spotSize(1)/voxSize(1));
    bmax(1) = max( min(bmax(1),imSize(1)), 1);
    bmax(2) = ceil(spotPosition(2) + spotSize(2)/voxSize(2));
    bmax(2) = max( min(bmax(2),imSize(2)), 1);
    if nDims == 3
        bmin(3) = ceil(spotPosition(3) - spotSize(3)/voxSize(3));
        bmin(3) = max( min(bmin(3),imSize(3)), 1);
        bmax(3) = ceil(spotPosition(3) + spotSize(3)/voxSize(3));
        bmax(3) = max( min(bmax(3),imSize(3)), 1);
    end
end