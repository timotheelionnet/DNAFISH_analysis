% builds an image (with arbitrary image_size) with randomly positioned 
% nspots spots, each the center of a 3D gaussian PSF (psfsigmaxy,psfsigmaz)

%% spots are always positioned at the center of the pixel they fall within.

%% parameters
img_size = [512,512,50]; %size of the stack: [nx ny nz] (integers).
nspots = 5; % number of spots in the image (integer)
psfsize = 21; %has to be odd number; the size of the kernel over which 
% the PSF is assumed to take non-negligible values
psfsigmaxy = 2; % sigma of gaussian PSF profile in the lateral plane
psfsigmaz = 5; % sigma of gaussian PSF profile in the axial direction

%% build artificial image and add center position of each spot
% zero evereywhere, except ones at the central pixel 
% of each randomly positioned spots
artificial_img = zeros(img_size);

spots_centers = rand(nspots,3);
spots_centers(:,1) = img_size(1)*spots_centers(:,1);
spots_centers(:,2) = img_size(2)*spots_centers(:,2);
spots_centers(:,3) = img_size(3)*spots_centers(:,3);
spots_centers = ceil(spots_centers);

%save position floats with loc3 file convention
save('testloc.loc3','spots_centers','-ascii');

%convert float positions with origin at image edge to pixel integers

spots_centers = sub2ind(img_size,...
    spots_centers(:,1),spots_centers(:,2),spots_centers(:,3));

artificial_img(spots_centers) = 1;

figure; imagesc( squeeze(max(artificial_img,[],3)) );
%you can also just generate an image in the workplace and name it
%artificial_img if you need deterministic spots, then run the rest of code
%section by section

%% convolve spot positions with PSF and save image
%build a grid of xyz coordinates
[y,x,z] = meshgrid( -floor(psfsize/2):floor(psfsize/2),...
    -floor(psfsize/2):floor(psfsize/2),...
    -floor(psfsize/2):floor(psfsize/2));

psf = exp( - (x.^2) / (2*psfsigmaxy^2)) ...
    .* exp( - (y.^2) / (2*psfsigmaxy^2)) ...
    .* exp( - (z.^2) / (2*psfsigmaz^2));

artificial_img = imfilter(artificial_img,psf,'conv');

figure; imagesc( squeeze(max(artificial_img,[],3)) );

save_as_tiff(artificial_img, 'testimg.tif',1);