function s = generate_gauss_spots_in_stack(pos,spotsize,imsize)
%All units in pixels
DefaultIntValue = 100;

ndims = numel(imsize);
s = zeros(imsize);
if ndims == 3
    if size(pos,2) <3
        disp('not enough columns');
    elseif size(pos,2) ==3
        pos(:,4) = DefaultIntValue;
    end   
elseif ndims ==2
    if size(pos,2) <2
        disp('not enough columns');
    elseif size(pos,2) ==2
        pos(:,3) = DefaultIntValue;
    end
end


for i=1:size(pos,1)
    
    %find coordinates of 3 sigma box around spot
    xmin = pos(i,1) - 3*spotsize(1);
    xmax = pos(i,1) + 3*spotsize(1);
    
    ymin = pos(i,2) - 3*spotsize(2);
    ymax = pos(i,2) + 3*spotsize(2);
    
    zmin = pos(i,3) - 3*spotsize(3);
    zmax = pos(i,3) + 3*spotsize(3);
    
    %make sure pixels are within the image
    xmin = round(max(1,xmin));
    xmin = round(min(imsize(1),xmin));
    xmax = round(max(1,xmax));
    xmax = round(min(imsize(1),xmax));
    
    ymin = round(max(1,ymin));
    ymin = round(min(imsize(2),ymin));
    ymax = round(max(1,ymax));
    ymax = round(min(imsize(2),ymax));
    
    zmin = round(max(1,zmin));
    zmin = round(min(imsize(3),zmin));
    zmax = round(max(1,zmax));
    zmax = round(min(imsize(3),zmax));
    
    [yy,xx,zz] = meshgrid(ymin:ymax,xmin:xmax,zmin:zmax);
    
    if min(size(xx))>0
        s(xmin:xmax,ymin:ymax,zmin:zmax) = s(xmin:xmax,ymin:ymax,zmin:zmax)+ ...
            pos(i,4)*...
            exp( - (xx-pos(i,1)).^2/(2*spotsize(1)^2))...
            .*exp( - (yy-pos(i,2)).^2/(2*spotsize(2)^2))...
            .*exp( - (zz-pos(i,3)).^2/(2*spotsize(3)^2));
    end
       
end


end

