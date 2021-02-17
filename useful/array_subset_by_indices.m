function out = array_subset_by_indices(in,index)

out = zeros(size(index,1),size(in,2));
for i =1:size(index,1)
    out(i,:) = in(index(i),:);    
end


end