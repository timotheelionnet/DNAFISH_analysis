function [i1f,i2f,d12,indices] = register_2images(i1,i2,fraction)

% registers two 2d-images i1 and i2. fraction is an integer that corresponds to the precision of the registration 
%(in fraction of pixel) 

%Outputs one subimage for each channel (i1f,i2f) that
% corresponds to the maximal overlap bewteen the both channels
%d12 is the translation offset between the two images in pixels

%% tests
if size(i1) ~= size(i2)
    disp('images have to be the same size');
    return
end

if ndims(i1) ~= 2
    disp('images should be 2D');
    return
end

%% finding out the translation offset
out12 = dftregistration(fft2(i1),fft2(i2),fraction);
%translation in x and y between two images
d12 = [out12(3),out12(4)];
clear('out12');

%% generating the registered overlapping images
if d12(1)>0
    i1f = i1(1+round(d12(1)): size(i1,1),:);
    indices(1,1,1) = 1+round(d12(1));
    indices(2,1,1) = size(i1,1);
    i2f = i2(1:size(i2,1)- round(d12(1)),:);
    indices(1,1,2) = 1; 
    indices(2,1,2) = size(i2,1)- round(d12(1));
else
    i1f = i1(1: size(i1,1)+round(d12(1)),:);
    indices(1,1,1) = 1;
    indices(2,1,1) = size(i1,1)+round(d12(1));
    i2f = i2(1- round(d12(1)):size(i2,1),:);
    indices(1,1,2) = 1- round(d12(1));
    indices(2,1,2) = size(i2,1);
end

if d12(2)>0
    i1f = i1f(:,1+round(d12(2)): size(i1f,2));
    indices(1,2,1) = 1+round(d12(2));
    indices(2,2,1) = size(i1f,2);
    i2f = i2f(:,1:size(i2f,2)- round(d12(2)));
    indices(1,2,2) = 1;
    indices(2,2,2) = size(i2f,2)- round(d12(2));
else
    i1f = i1f(:,1: size(i1f,2)+round(d12(2)));
    indices(1,2,1) = 1;
    indices(2,2,1) = size(i1f,2)+round(d12(2));
    i2f = i2f(:,1- round(d12(2)):size(i2f,2));
    indices(1,2,2) = 1- round(d12(2));
    indices(2,2,2) = size(i2f,2);
end

end
