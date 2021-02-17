function params = quadratic_estimate3D_old(spots,dx)
%fits dx = ax2+by2+cxy+dx+ey+f
%spots = [x,y,z]

%map coordinates data to [0,1] to avoid diverging sums:
minx = min(spots(:,1));
miny = min(spots(:,2));
%minz = min(spots(:,3));

xnorm = max(spots(:,1))-minx;
ynorm = max(spots(:,2))-miny;
%znorm = max(spots(:,3))-minz;
spots(:,1) = (spots(:,1)-minx)/xnorm;
spots(:,2) = (spots(:,2)-miny)/ynorm;
%spots(:,3) = (spots(:,3)-minz)/znorm;

%compute sums:
Sx = sum(spots(:,1));
Sy = sum(spots(:,2));
%Sz = sum(spots(:,3));
Sw = sum(dx);

Sx2 = sum(spots(:,1).^2);
Sy2 = sum(spots(:,2).^2);
%Sz2 = sum(spots(:,3).^2);
Sxy = sum(spots(:,1).*spots(:,2));
%Sxz = sum(spots(:,1).*spots(:,3));
%Syz = sum(spots(:,2).*spots(:,3));
Sxw = sum(spots(:,1).*dx);
Syw = sum(spots(:,2).*dx);
%Szw = sum(spots(:,3).*dx(:,3));

Sx2w = sum(spots(:,1).^2.*dx);
Sy2w = sum(spots(:,2).^2.*dx);
Sxyw = sum(spots(:,1).*spots(:,2).*dx);
%Sz2w = sum(spots(:,3).^2.*dx(:,3));

Sx3 = sum(spots(:,1).^3);
Sy3 = sum(spots(:,2).^3);
%Sz3 = sum(spots(:,3).^3);
Sx2y = sum(spots(:,1).^2.*spots(:,2));
%Sx2z = sum(spots(:,1).^2.*spots(:,3));
%Sy2z = sum(spots(:,2).^2.*spots(:,3));
Sxy2 = sum(spots(:,1).*spots(:,2).^2);
%Sxz2 = sum(spots(:,1).*spots(:,3).^2);
%Syz2 = sum(spots(:,2).*spots(:,3).^2);

Sx4 = sum(spots(:,1).^4);
Sy4 = sum(spots(:,2).^4);
%Sz4 = sum(spots(:,3).^4);
Sx2y2 = sum(spots(:,1).^2.*spots(:,2).^2);
%Sx2z2 = sum(spots(:,1).^2.*spots(:,3).^2);
%Sy2z2 = sum(spots(:,2).^2.*spots(:,3).^2);
Sx3y = sum(spots(:,1).^3.*spots(:,2));
Sxy3 = sum(spots(:,1).*spots(:,2).^3);

N = size(spots,1);

%6x6 matrix to invert
M =[Sx4,  Sx2y2,Sx3y,Sx3, Sx2y,Sx2;...
    Sx2y2,Sy4,  Sxy3,Sxy2,Sy3, Sy2;...
    Sx3y,Sxy3,Sx2y2,  Sx2y,Sxy2,Sxy;...
    Sx3  ,Sxy2 ,Sx2y, Sx2, Sxy, Sx;...
    Sx2y, Sy3,  Sxy2, Sxy, Sy2, Sy;...
    Sx2,  Sy2,  Sxy,  Sx,  Sy,  N];

p =  M^(-1)*[Sx2w;Sy2w;Sxyw;Sxw;Syw;Sw]; 

%correcting parameters values for normalization
params = zeros(6,1);
params(1) = p(1)/xnorm^2;
params(2) = p(2)/ynorm^2;
params(3) = p(3)/(xnorm*ynorm);

params(4) = p(4)/xnorm -2*p(1)*minx/xnorm^2 - p(3)*miny/(xnorm*ynorm);
params(5) = p(5)/ynorm -2*p(2)*miny/ynorm^2 - p(3)*minx/(xnorm*ynorm);

params(6) = p(6)+ p(1)*minx^2/xnorm^2 - p(4)*minx/xnorm...
            + p(2)*miny^2/ynorm^2 - p(5)*miny/ynorm...
            + p(3)*minx*miny/(xnorm*ynorm);

end

