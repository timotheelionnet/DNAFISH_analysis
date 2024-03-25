function dist = get_dist_matrix(pos1,pos2)
    %pos1 = [x y z] (m rows)        pos2 = [x y z] (n rows)
    %dist is the mxn matrix of the distances between spots.

    m = size(pos1,1); n = size(pos2,1);
    
    dx = repmat(pos1(:,1),1,n) - repmat(pos2(:,1)',m,1);
    dx = dx.^2;
    dy = repmat(pos1(:,2),1,n) - repmat(pos2(:,2)',m,1);
    dy = dy.^2;
    dz = repmat(pos1(:,3),1,n) - repmat(pos2(:,3)',m,1);
    dz = dz.^2;
    
    
    dist = dx + dy + dz;
    dist = sqrt(dist);
    clear('dx','dy','dz','m','n');
end