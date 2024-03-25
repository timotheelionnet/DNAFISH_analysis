function stackx = convert_coordinates_to_stack_squares(varargin)

%takes the positions input as a list and builds a stackx (of dimensions given as input) 
%where hollow squares mark each position

%convert_coordinates_to_stackx([x y z I],nx,ny,nz)
%convert_coordinates_to_stackx([x y z],nx,ny,nz)
%convert_coordinates_to_stackx([x y z I],nx,ny,nz,valmode)

%convert_coordinates_to_stackx([x y I],nx,ny)
%convert_coordinates_to_stackx([x y],nx,ny)
%convert_coordinates_to_stackx([x y I],nx,ny,valmode)

%valmode optional argument
%valmode = 'value' (default): I build a stackx with the values corresponding to the last
%column of the array [x y z I] or [x y I]

%valmode = 'ones': I put ones in the stackx at the positions indicated by the
%input array [x y z I] or [x y I]

square_size = 3; % the radius of the square around each pixel


if nargin == 0
    stackx = 0;
    return;
end

if nargin >=1
    coor = varargin{1};
    %% 3d case
    if nargin == 5 || (nargin == 4 && ~ischar(varargin{4}) )  %3d case
        nx = varargin{2}; ny = varargin{3}; nz = varargin{4};
        if nargin ==5
            valmode = varargin{5};
        else
            if size(coor,2) >= 4
                valmode = 'value';
            elseif size(coor,2) == 3
                valmode = 'ones';
            else
                disp(['3d convert_coordinates_to_stackx: input array has the wrong column number: ' num2str(size(coor,2))]);
            end
        end

        stackx = fill_stack_3d_squares(coor,nx,ny,nz,valmode,square_size);
        return;

    %% 2d case
    elseif nargin == 3 || (nargin == 4 && ischar(varargin{4}) )  
        nx = varargin{2}; ny = varargin{3};
        if nargin ==4
            valmode = varargin{4};
        else

            if size(coor,2) >= 3
                valmode = 'value';
            elseif size(coor,2) == 2
                valmode = 'ones';
            else
                disp(['2d convert_coordinates_to_stackx: input array has the wrong column number: ' num2str(size(coor,2))]);
            end
        end

        stackx = fill_stack_2d_squares(coor,nx,ny,valmode,square_size);
        return

    %% wrong number/format of arguments    
    else        
        return
    end
end
clear('valmode','nx','ny','nz');
end


function stackx = fill_stack_3d_squares(coor,nx,ny,nz,valmode,square_size)
% fill in a 3D stackx (1:nx,1:ny,1:nz) with zeros except at positions indicated by coor = [x y z I] 
%or [x y z].
%if input positions falls outside of the size of the array, I change the coordinate
%to zero or nx/ny/nz
    
    if size(coor,2) == 3
        stackx = zeros(nx,ny,nz,'uint8');
    elseif size(coor,2) == 4
        stackx = zeros(nx,ny,nz,class(coor));
    end

    if ndims(coor)<2
        return
    end

    %I impose that coordinates fall within the destination array
    coor(:,1) = max(1,min(uint16(ceil(coor(:,1))),nx));
    coor(:,2) = max(1,min(uint16(ceil(coor(:,2))),ny));
    coor(:,3) = max(1,min(uint16(ceil(coor(:,3))),nz));

    for i = 1:size(coor,1)
        if size(coor,2)>=4 && strcmp(valmode,'value') 
            stackx = draw_square_around_pixel_3d(stackx,coor(i,1),coor(i,2),coor(i,3),nx,ny,nz,coor(i,4),square_size);    
        elseif strcmp(valmode,'ones')  
            stackx = draw_square_around_pixel_3d(stackx,coor(i,1),coor(i,2),coor(i,3),nx,ny,nz,1,square_size);   
        else
            stackx = draw_square_around_pixel_3d(stackx,coor(i,1),coor(i,2),coor(i,3),nx,ny,nz,1,square_size);
            
        end
    end
    clear('coor','i');
end


function stackx = draw_square_around_pixel_3d(stackx,xp,yp,zp,nx,ny,nz,value,square_size)

if zp>0 && zp<=nz
    
    for i = xp-square_size:xp+square_size
        j1 = yp - square_size;
        if i>0 && i<=nx && j1>0 && j1<ny 
            stackx(i,j1,zp) = value;
        end 
        j2 = yp + square_size;
        if i>0 && i<=nx && j2>0 && j2<ny 
            stackx(i,j2,zp) = value;
        end 
    end
    
    for j = yp-square_size+1:yp+square_size-1
        i1 = xp - square_size;
        if j>0 && j<=ny && i1>0 && i1<nx
            stackx(i1,j,zp) = value;
        end 
        i2 = xp + square_size;
        if j>0 && j<=ny && i2>0 && i2<nx 
            stackx(i2,j,zp) = value;
        end 
    end   
end

clear('i','i1','i2','j','j1','j2');
end



function stackx = fill_stack_2d_squares(coor,nx,ny,valmode,square_size)
% fill in a stackx [nx ny nz] with zeros except at positions indicated by coor = [x y z I] 
%or [x y z].
%if input positions falls outside of the size of the array, I change the coordinate
%to zero or nx/ny/nz
    if size(coor,2) == 2
        stackx = zeros(nx,ny,'uint8');
    elseif size(coor,2) == 3
        stackx = zeros(nx,ny,class(coor));
    end

    if ndims(coor)<2
        return
    end  

    coor(:,1) = max(1,min(uint16(ceil(coor(:,1))),nx));
    coor(:,2) = max(1,min(uint16(ceil(coor(:,2))),ny));
    
    for i = 1:size(coor,1)
        if size(coor,2)>=3 && strcmp(valmode,'value')            
            stackx = draw_square_around_pixel_2d(stackx,coor(i,1),coor(i,2),nx,ny,coor(i,3),square_size);    
        elseif strcmp(valmode,'ones')  
            stackx = draw_square_around_pixel_2d(stackx,coor(i,1),coor(i,2),nx,ny,1,square_size);  
        else
            stackx = draw_square_around_pixel_2d(stackx,coor(i,1),coor(i,2),nx,ny,1,square_size);  
        end
    end
    clear('i','coor');
end

function stackx = draw_square_around_pixel_2d(stackx,xp,yp,nx,ny,value,square_size)

for i = xp-square_size:xp+square_size
        j1 = yp - square_size;
        if i>0 && i<=nx && j1>0 && j1<ny 
            stackx(i,j1) = value;
        end 
        j2 = yp + square_size;
        if i>0 && i<=nx && j2>0 && j2<ny 
            stackx(i,j2) = value;
        end 
end
    
for j = yp-square_size+1:yp+square_size-1
        i1 = xp - square_size;
        if j>0 && j<=ny && i1>0 && i1<nx
            stackx(i1,j) = value;
        end 
        i2 = xp + square_size;
        if j>0 && j<=ny && i2>0 && i2<nx 
            stackx(i2,j) = value;
        end 
end   


clear('i','i1','i2','j','j1','j2');
end
