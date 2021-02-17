%%
%load stacks in each channel
s1 = tiffread5('/Users/lionnett/Documents/Salina Long/20160613_mi4_data/C1-mi4Halo_Gad1ax488(green)_vGlutCy3(red)_HaloJF646(blue)-fullstack.tif'); %blue - HALO
s2 = tiffread5('/Users/lionnett/Documents/Salina Long/20160613_mi4_data/C2-mi4Halo_Gad1ax488(green)_vGlutCy3(red)_HaloJF646(blue)-fullstack.tif'); %red - vGlut
s3 = tiffread5('/Users/lionnett/Documents/Salina Long/20160613_mi4_data/C3-mi4Halo_Gad1ax488(green)_vGlutCy3(red)_HaloJF646(blue)-fullstack.tif'); %green - Gad1

%compute mean and SD for each channel
mu1 = mean(double(s1(:))); 
sig1 = std(double(s1(:)));

mu2 = mean(double(s2(:))); 
sig2 = std(double(s2(:)));

mu3 = mean(double(s3(:))); 
sig3 = std(double(s3(:)));

%% compute threshold
n = 3; %number of STDs away from mean image level

%indices of background pixels
idx = double( (s1(:) < mu1 + n*sig1) ) .* double( (s2(:) < mu2 + n*sig2) ) .* double( (s3(:) < mu3 + n*sig3) );

%convert to logical, invert to get the non-background pixels
idx = ~logical(idx);

%collect all non background pixels
s1f = double(s1(idx));
s2f = double(s2(idx));
s3f = double(s3(idx));

%map intensity distributions to mean = 0, mean + STD = 1
s1fnorm = (s1f - mu1)/sig1;
s2fnorm = (s2f - mu2)/sig2;
s3fnorm = (s3f - mu3)/sig3;

%% plot as spot cloud (not great because too many spots cause matlab too become too slow)
plot_spots = 0; %to enable plotting all the spots, set variable to 1; to disable, set to zero

if plot_spots
    figure('Name','1 2'); 
    plot(s1fnorm(:),s2fnorm(:)); 

    figure('Name','1 3'); 
    plot(s1fnorm(:),s3fnorm(:)); 

    figure('Name','2 3'); 
    plot(s2fnorm(:),s3fnorm(:));
end

%% compute 2D histograms (better representation)
%number of bins in the 2D histogram
nbins = [20,20];

[h12,c12] = hist3([s1fnorm(:),s2fnorm(:)],'Edges',{0:2:50,0:2:50});
figure('Name','1 2'); 
imagesc(log(h12+1));
axis xy

[h13,c13] = hist3([s1fnorm(:),s3fnorm(:)],'Edges',{0:2:50,0:2:50});
figure('Name','1 3'); 
imagesc(log(h13+1));
axis xy

[h23,c23] = hist3([s2fnorm(:),s3fnorm(:)],'Edges',{0:2:50,0:2:50});
figure('Name','2 3'); 
imagesc(log(h23+1));
axis xy

% figure('Name','1 3'); 
% plot(s1fnorm(:),s3fnorm(:)); 
% 
% figure('Name','2 3'); 
% plot(s2fnorm(:),s3fnorm(:));


