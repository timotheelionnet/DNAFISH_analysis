function coor = convert_stack_to_coordinates(stack)

coor = zeros(1000,4);

%% 3d case
if ndims(stack) == 3    
    if sum(sum(sum(stack)))==0
        disp('no nonzero points in your stack');
        return
    else
        tmp = sum(sum(sum(stack)));
    end

    i = 1; tsum = 0; n=1;
    while i <= size(stack,1) && tsum < tmp
        if sum(sum(stack(i,:,:))) ~=0
            tmpi = sum(sum(stack(i,:,:)));
            j=1; sumi = 0;
            while j <= size(stack,3) && sumi < tmpi
                if sum(stack(i,:,j)) ~=0
                    tmpj = sum(stack(i,:,j));
                    k=1; sumj = 0;
                    while k<=size(stack,2) && sumj < tmpj
                        if stack(i,k,j) ~=0
                            coor(n,1) = i; coor(n,2) = k; coor(n,3) = j; coor(n,4) = stack(i,k,j);
                            sumj = sumj + stack(i,k,j);
                            sumi = sumi + stack(i,k,j);
                            tsum = tsum + stack(i,k,j);
                            n = n+1;
                        end
                        k = k+1;
                    end
                end
                j=j+1;
            end
        else 
        end
        i = i+1;
    end

    coor = coor(1:n-1,:);
    if sum(coor(:,4)) ~= tsum
        disp('there was a problem in the conversion of the stack to coordinates');
    end

%% 2d case    
elseif ndims(stack) == 2
    if sum(sum(stack)) ==0
        disp('no nonzero points in your stack');
        return
    else
        tmp = sum(sum(stack));
    end

    i = 1; tsum = 0; n=1;
    while i <= size(stack,1) && tsum < tmp
        if sum(stack(i,:)) ~=0
            tmpi = sum(stack(i,:));
            j=1; sumi = 0;
            while j <= size(stack,2) && sumi < tmpi
                if stack(i,j) ~=0
                    coor(n,1) = i; coor(n,2) = j; coor(n,3) = stack(i,j);
                    sumi = sumi + stack(i,j);
                    tsum = tsum + stack(i,j);
                    n = n+1;
                end
                j = j+1;
            end
        end
        i = i+1;
    end
    coor = coor(1:n-1,:);
    if sum(coor(:,3)) ~= tsum
        disp('there was a problem in the conversion of the stack to coordinates');
    end
end
end
