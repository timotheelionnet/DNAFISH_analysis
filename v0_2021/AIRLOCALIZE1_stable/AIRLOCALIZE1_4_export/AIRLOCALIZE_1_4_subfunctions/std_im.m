function sd = std_im(img)
    %2d std
    [nx ny] = size(img);
    sd = std(double(reshape(img,1,nx*ny)));
    clear('nx','ny');
end