function out = subsample_data_set(arr,npts)

[nx,ny] = size(arr);
nout = floor(nx/npts);
if nout < 1
    disp(['your array size = ' num2str(nx) ' is smaller than your sampling interval ' num2str(npts)]);
    out = arr;
    return;
end
out = zeros(nout,ny);

out(1:nout,1:ny) = arr(1:npts:(nout-1)*npts+1,1:ny);

end