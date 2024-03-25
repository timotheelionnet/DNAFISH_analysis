%% import the spots
%col 1: x in pixels
%col 2: y in pixels
%col 3: z in pixels
%col 4: Intensity 
%col 5: Intensity residuals (not used)
%col 6: frame (not used)
p1 = r1.final_pix;
p2 = r2.final_pix;

%% remove artefactual spots 
%(negative intensities can happen when spot detection didnt converge)
p1 = p1(p1(:,4)>0,:);
p2 = p2(p2(:,4)>0,:);

%% plot Intensity Distribution
nbins_I = 1000;
I1max = p1(1,4);
I1_edges = 0:I1max/nbins_I:I1max;
[nI1,I1_edges] = histcounts(p1(:,4),I1_edges);

I2max = p2(1,4);
I2_edges = 0:I2max/nbins_I:I2max;
[nI2,I2_edges] = histcounts(p2(:,4),I2_edges);

figure; hold;
plot((I1_edges(2:end) + I1_edges(1:end-1))/2,nI1);
plot((I2_edges(2:end) + I2_edges(1:end-1))/2,nI2);

%% compute distance matrices
n1 = size(p1,1);
n2 = size(p2,1);

%between the two channels
dx12 = repmat(p1(:,1),1,n2) - repmat(p2(:,1)',n1,1);
dy12 = repmat(p1(:,2),1,n2) - repmat(p2(:,2)',n1,1);
dz12 = repmat(p1(:,3),1,n2) - repmat(p2(:,3)',n1,1);
dr12 = sqrt(dx12.^2 + dy12.^2 + dz12.^2);

%channel 1 alone
dx11 = repmat(p1(:,1),1,n1) - repmat(p1(:,1)',n1,1);
dy11 = repmat(p1(:,2),1,n1) - repmat(p1(:,2)',n1,1);
dz11 = repmat(p1(:,3),1,n1) - repmat(p1(:,3)',n1,1);
dr11 = sqrt(dx11.^2 + dy11.^2 + dz11.^2);
%keeping only the non redundant pairs
dr11 = dr11(find(~tril(ones(size(dr11)))));

dx22 = repmat(p2(:,1),1,n2) - repmat(p2(:,1)',n2,1);
dy22 = repmat(p2(:,2),1,n2) - repmat(p2(:,2)',n2,1);
dz22 = repmat(p2(:,3),1,n2) - repmat(p2(:,3)',n2,1);
dr22 = sqrt(dx22.^2 + dy22.^2 + dz22.^2);
%keeping only the non redundant pairs
dr22 = dr22(find(~tril(ones(size(dr22)))));

%% plot absolute offset in each dimensions
x_edges = -100:0.1:100;
[nx,x_edges] = histcounts(dx12(:),x_edges);
[ny,y_edges] = histcounts(dy12(:),x_edges);
[nz,z_edges] = histcounts(dz12(:),x_edges);

figure; hold;
plot((x_edges(2:end) + x_edges(1:end-1))/2,nx);
plot((y_edges(2:end) + y_edges(1:end-1))/2,ny);
plot((z_edges(2:end) + z_edges(1:end-1))/2,nz);

%% plot the distances between spots
r_edges = 0:0.02:100;
[nr12,r_edges] = histcounts(dr12(:),r_edges);
[nr11,r_edges] = histcounts(dr11(:),r_edges);
[nr22,r_edges] = histcounts(dr22(:),r_edges);
vol = 4*pi/3*r_edges(2:end).^3 - 4*pi/3*r_edges(1:end-1).^3;
figure; hold;
plot((r_edges(2:end) + r_edges(1:end-1))/2,nr12./vol);
plot((r_edges(2:end) + r_edges(1:end-1))/2,nr11./vol);
plot((r_edges(2:end) + r_edges(1:end-1))/2,nr22./vol);

%% pair spots with colocalizing spot
dist_threshold = 3;
Idx=zeros(n1,3);
for i=1:size(dr12,1)
    [Idx(i,2),Idx(i,1)] = min(dr(i,:));
    if(Idx(i,2)>dist_threshold)
        Idx(i,3) = 0;
    else
        Idx(i,3) = Idx(i,1);
    end
end


    
    
    