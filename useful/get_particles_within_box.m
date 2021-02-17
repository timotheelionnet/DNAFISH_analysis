function opts = get_particles_within_box(ipts,xmin,xmax,ymin,ymax,zmin,zmax)

indices = (ipts(:,1)<xmax ).*(ipts(:,1)>xmin ).*(ipts(:,2)<ymax ).*(ipts(:,2)>ymin ).*(ipts(:,3)<zmax ).*(ipts(:,3)>zmin );

opts = ipts(find(indices),:);
clear('indices');
end