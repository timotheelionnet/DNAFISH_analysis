function st2 = translate_stack(st1,dx,dy,dz)
%st has same size as st1, the points in st2 are obtained from st1 by
%translation of dx,dy,dz 
%padding is performed with the median of st1 as default value.
[nx ny nz] = size(st1);
bg = median_stack(st1);
dx = round(dx); dy = round(dy); dz = round(dz);
st2 = bg*ones(nx,ny,nz);

if(dx>=0)
    x1min = 1; x1max = nx-dx; x2min = 1+dx; x2max = nx;
else
    x1min = 1-dx; x1max = nx; x2min = 1; x2max = nx + dx;
end
if(dy>=0)
    y1min = 1; y1max = ny-dy; y2min = 1+dy; y2max = ny;
else
    y1min = 1-dy; y1max = ny; y2min = 1; y2max = ny + dy;
end
if(dz>=0)
    z1min = 1; z1max = nz-dz; z2min = 1+dz; z2max = nz;
else
    z1min = 1-dz; z1max = nz; z2min = 1; z2max = nz + dz;
end
st2(x2min:x2max,y2min:y2max,z2min:z2max) = st1(x1min:x1max,y1min:y1max,z1min:z1max);
end