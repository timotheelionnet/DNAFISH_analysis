function stack2 = generate_overlay_channel(varargin)
    
    args = varargin{1};
    nargs = size(varargin,2);
    [nx ny nz] = size(args{1});
    
    if nargs>=2
        stack2 = args{2};
        if ndims(stack2) ~= 3 || size(stack2,3)==1
            fprintf('I like my arguments to be 3 dimensional, babe.\n');
            return;
        end
    else    stack2 = zeros(nx,ny,nz);
    end
    
    if [nx ny nz] ~= size(stack2)
        fprintf('I like my stacks to be the same size, babe.\n');
        return;
    end
    
    if nargs >= 3  
        if strcmp('small',args{3})
        elseif strcmp('large',args{3})
            stack2 = dilate(stack2,'cube');
        elseif strcmp('mix',args{3})
            stack2 = dilate(stack2,'mix');
        elseif strcmp('square',args{3})
            stack2 = dilate(stack2,'square');
        end 
    else
        if nargs>=2
        stack2 = dilate(stack2,'square');
        end
    end
    
    stack2(:,:,:) = max(  min( stack2(:,:,:) , 1 ) , 0) ; %overlay data should be between 0 and 1. Convolution can lead to >1 values.
    clear('args','nargs','nx','ny','nz');
end