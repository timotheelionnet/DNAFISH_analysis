function dist = get_dist_matrix2d(pos1,pos2)
    %pos1 = [x y] (m rows)        pos2 = [x y] (n rows)
    %dist is the mxn matrix of the distances between spots.

    m = size(pos1,1); n = size(pos2,1);
    
    dx = repmat(pos1(:,1),1,n) - repmat(pos2(:,1)',m,1);
    dx = dx.^2;
    dy = repmat(pos1(:,2),1,n) - repmat(pos2(:,2)',m,1);
    dy = dy.^2;
    dist = dx + dy;
    dist = sqrt(dist);
    clear('dx','dy','dz','m','n');
end