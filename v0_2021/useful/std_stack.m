function sd = std_stack(arr)
    %does not use the buit in matlab function that has memory leaks
    [nx ny nz] = size(arr);
    m = zeros(1,nz);
    m2 = zeros(1,nz);
    
    for i =1:nz
        m(1,i) = mean(reshape(arr(:,:,i),nx*ny,1));
    end
    m = mean(m,2);
    for i=1:nz
        m2(1,i) = mean(reshape((arr(:,:,i) - m).^2,nx*ny,1));
    end
    sd = mean(m2,2);
    
    clear('nx','ny','nz','m','m2','i');
end