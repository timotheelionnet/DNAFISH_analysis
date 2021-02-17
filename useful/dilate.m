function dil = dilate(stack,shape)

if ndims(stack) ==3
    [nx ny nz] = size(stack);
    dil = zeros(nx,ny,nz);
elseif ndims(stack) == 2
    [nx ny] = size(stack);
    dil = zeros(nx,ny);
else return
end


 if strcmp(shape,'cross')
     B = generate_cross_matrix(ndims(stack));
     dil = convn(stack,B,'same');
     
 elseif strcmp(shape,'cube')   
     if ndims(stack) == 3
     B = ones(3,3,3);
     elseif ndims(stack) == 2
         B = ones(3,3);
     end
     dil = convn(stack,B,'same');
     
 elseif strcmp(shape,'mix')
     B = generate_mix_matrix(nd);
     dil = convn(stack,B,'same');
     
 elseif strcmp(shape,'square')   
     if ndims(stack) == 3
         B = ones(3,3,1);
     elseif ndims(stack) == 2
         B = ones(3,3);
     end
     dil = convn(stack,B,'same');
 end

dil(:,:,:) = min(dil(:,:,:),1);
dil(:,:,:) = max(dil(:,:,:),0);

clear('stack','B');
end

function B = generate_cross_matrix(nd)
     if (nd == 3)
         B = zeros(3,3,3);
         B(2,2,2) = 1;
         B(2,2,1) = 1;
         B(2,2,3) = 1;
         B(2,1,2) = 1;
         B(2,3,2) = 1;
         B(1,2,2) = 1;
         B(3,2,2) = 1;
     elseif nd == 2
         B = zeros(3,3);
         B(2,2) = 1;
         B(2,1) = 1;
         B(2,3) = 1;
         B(1,2) = 1;
         B(3,2) = 1;
     end
end

function B = generate_mix_matrix(nd)
     if(nd == 3)
         B = zeros(3,3,3);
         B(:,:,1) = [0,1,0;1,1,1;0,1,0];
         B(:,:,2) = ones(3,3);
         B(:,:,3) = [1,0,1;0,1,0;1,0,1];
         
     elseif (nd ==2)
         B = zeros(3,3);
         B(2,2) = 1;
         B(2,1) = 1;
         B(2,3) = 1;
         B(1,2) = 1;
         B(3,2) = 1;
     end
end