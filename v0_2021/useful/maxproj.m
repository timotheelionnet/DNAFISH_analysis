function proj = maxproj(pic)
%makes the max projection of an image input as a 3D stack.
proj(:,:) = double(squeeze(max(pic,[],3)));
end
