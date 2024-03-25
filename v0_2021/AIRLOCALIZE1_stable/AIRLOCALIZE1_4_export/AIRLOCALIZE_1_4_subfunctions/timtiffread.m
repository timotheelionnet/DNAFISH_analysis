function [stack, filename] = timtiffread(varargin)

if nargin <1
    %disp('before tiffread');
    %memory
    [stack, nbImages, filename] = tiffread5;    
else
    filename = varargin{1};
    [stack, nbImages, filename] = tiffread5(filename);
end

% %% the following section is obsolete with tiffread5 
% %because tiffread5 outputs a 3d array, not a structure. The
% %reason is that I couldnt get the following part to operate without
% %fragmenting the memory too much.
% if(isstruct(stack))
%     stack = uint16(stack);
%     nx = size(stack(1).data,1);
%     ny = size(stack(1).data,2);
% 
%     stack_arr = zeros(nx,ny,nbImages,'uint16');
% 
%     for i=1:nbImages
%         stack_arr(:,:,i) = stack(i).data;
%     end
%     stack = stack_arr;    
% end

%%
clear('stack_arr','nbImages','nx','ny','i');
%disp('after tiffread');
%memory
end

