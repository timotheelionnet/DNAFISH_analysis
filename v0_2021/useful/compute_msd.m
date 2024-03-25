function [msd,dt,n] = compute_msd(X)

n_detections = size(X, 1);
dt = 1:n_detections - 1;
msd = zeros(numel(dt));
n = zeros(numel(dt));
for j = 1 : numel(dt)   
    dr = []; 
    k = 0;
    tmp_msd = NaN;
    for i=1:n_detections - dt(j)
        dr2 = sum(X(i,:) - X(i+dt(j),:),2).^2;
        if ~isnan(dr)
            k = k+1;
            if k == 1;
                tmp_msd = dr2;
            else
                tmp_msd = tmp_msd + dr2;
            end
        end
    end
    msd(j) = tmp_msd/k;
    n(j) = k;
e