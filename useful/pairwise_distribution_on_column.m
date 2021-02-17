function pdis = pairwise_distribution_on_column(arr,colnum)

npts = size(arr,1);

arr1 = repmat(arr(:,colnum),1,npts);
arr2 = repmat(arr(:,colnum)',npts,1);

diff = abs(arr1 - arr2);
diff = tril(diff,-1);
pdis = nonzeros(diff);
clear('arr1','arr2','diff','npts','arr','colnum');
end