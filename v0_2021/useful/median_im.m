function med = median_im(img)
    %2d median
    [nx ny] = size(img);
    med = median(double(reshape(img,1,nx*ny)));
    clear('nx','ny');
end