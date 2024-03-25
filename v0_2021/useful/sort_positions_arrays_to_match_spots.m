function list2sorted = sort_positions_arrays_to_match_spots(varargin)
%inputs two lists of spots formatted as follows [x,y,z] or [x,y]:
%sort_positions_arrays_to_match_spots(list1,list2)
%a third optional argument can be input to force the dimensionality, in
%case the arrays have more columns than their number of dimensions

list1 = varargin{1};
list2 = varargin{2};
ndims = 0;
if nargin >=3
    ndims = varargin{3};
else
    ndims = min(3,size(list1,2));
end


if size(list2,1) == 1 || ndims == 0;
    list2sorted = list2;
    return
end
%% distance matrix
[xx1,xx2] = meshgrid(list2(:,1),list1(:,1));
xx1 = xx1-xx2; 
clear('xx2');
xx1 = xx1.*xx1;


[yy1,yy2] = meshgrid(list2(:,2),list1(:,2));
yy1 = yy1-yy2; 
clear('yy2');
yy1 = yy1.*yy1;

if ndims == 3
    [zz1,zz2] = meshgrid(list2(:,3),list1(:,3));
    zz1 = zz1-zz2; 
    clear('zz2');
    zz1 = zz1.*zz1;
    dist = zz1+xx1+yy1;
else
    dist = xx1+yy1;
end

clear('zz1','xx1','yy1');

[min12,imin12] = min(dist);
[min21,imin21] = min(dist,[],2);
clear('dist');

%% matching nearest neighbours
list2sorted = zeros(size(list2,1),size(list2,2)+1);
j=1;

if size(list1,1) == 1
    for i = 1:size(list2,1)
        if i == imin21(1)
            list2sorted(1,1:size(list2,2))= list2(i,1:size(list2,2));
            list2sorted(1,size(list2,2)+1) = i;
        else
            list2sorted(size(list1,1)+j,1:size(list2,2))= list2(i,1:size(list2,2));
            list2sorted(size(list1,1)+j,size(list2,2)+1) = i;
            j=j+1;
        end
    end
    
else

    for i = 1:size(list2,1)
        if i == imin21(imin12(i))
            list2sorted(imin12(i),1:size(list2,2))= list2(i,1:size(list2,2));
            list2sorted(imin12(i),size(list2,2)+1) = i;
        else
            list2sorted(size(list1,1)+j,1:size(list2,2))= list2(i,1:size(list2,2));
            list2sorted(size(list1,1)+j,size(list2,2)+1) = i;
            j=j+1;
        end
    end
end
%% housekeeping
clear('x1','y1','z1','x2','y2','z2','xx1','xx2','yy1','yy2','zz1','zz2','dx','dy','dz',...
    'i','j','min12','imin12','min21','imin21','dist','ndims');
end