function r2 = Translate_rotate_dataset(r1,r0,thetax,thetay,thetaz)
%r1 = data, N x 3 array
%with [x,y,z ] columns

%r0 = [x0,y0,z0] translation vector

%thetax, thetay, thetaz Euler angles

Rx =   [1,              0,              0;...
        0,              cos(thetax),    -sin(thetax);...
        0,              sin(thetax),     cos(thetax)];

Ry =   [cos(thetay),    0,              sin(thetay);...
        0,              1,              0;...
        -sin(thetay),   0,              cos(thetay)];

Rz =   [cos(thetaz),    -sin(thetaz),   0;...
        sin(thetaz),    cos(thetaz),    0;...
        0,              0,              1];

r0 = repmat(r0,size(r1,1),1); 
    
r1 = Rz*Ry*Rx*r1';
r2 = (r1 + r0')';





end

