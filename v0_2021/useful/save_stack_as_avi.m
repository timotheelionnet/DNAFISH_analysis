function save_stack_as_avi(stack,fname,normalize)
    %stack must be a z-stack (must be 3D)
    %converts into 8 bit by scaling the dynamic range of the stack onto the
    %[0,255] interval
    %saves as avi
    
    %choose normalization between by slice and global over the entire stack
    %normalize = 'global'
    %or
    %normalize = 'slice'
    
    %normalize = 'global';
    %normalize = 'slice';
    
    switch normalize
        case 'global'
            stack = double(stack);
            mov = immovie(((reshape(stack,size(stack,1),size(stack,2),1,size(stack,3)) - min(stack(:)))*255)/(max(stack(:)) - min(stack(:)))+1,colormap(gray(256)));
            
        case 'slice'
            minz = min(min(stack,[],2),[],1);
            normz = max(max(stack,[],2),[],1) - minz;

            minz = reshape(minz,1,1,numel(minz));
            normz = reshape(normz,1,1,numel(normz));

            stack = stack - repmat(minz,size(stack,1),size(stack,2),1);
            normz = repmat(normz,size(stack,1),size(stack,2),1);
            
            stack = double(stack)*255./double(normz) +1;
            mov = immovie(reshape(stack,size(stack,1),size(stack,2),1,size(stack,3)),colormap(gray(256)));
    end
    
    YourVideo = VideoWriter(fname); 
    
    open(YourVideo); 
    
    writeVideo(YourVideo,mov);
    
    close(YourVideo);
    
    
end