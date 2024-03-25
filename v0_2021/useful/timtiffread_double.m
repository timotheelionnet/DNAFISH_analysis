function [stack_arr, filename] = timtiffread_double(varargin)

if nargin <1
    [stack, nbImages, filename] = tiffread2;
else
    [stack, nbImages, filename] = tiffread2(varargin{1});
end

nx = size(stack(1).data,1);
ny = size(stack(1).data,2);

stack_arr = zeros(nx,ny,nbImages);

for i=1:nbImages
stack_arr(:,:,i) = stack(i).data;
end

clear('stack','nbImages','nx','ny');

end