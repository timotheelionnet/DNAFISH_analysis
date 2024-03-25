function proj = maxproj_struc(pic)
    %makes the max projection of an image input as a structure pic.

    p =pic.ndim;

    if pic.ndim == 2
         disp('only one z plane in your image');
         proj = pic.data;

    elseif pic.ndim == 3
        if pic.nimg == 1
            im = pic.data;
        else
            im = pic.data{1};
        end
        proj = squeeze(max(im,[],3));
    else
        proj = 0;
    end
    
    
end
