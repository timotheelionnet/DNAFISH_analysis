function out = file_exists(fname,sourcedir)
    list = dir(sourcedir);
    if max(strcmp(fname,{list.name}')) ~= 0
        out = 1;
    else
        out = 0;
    end
    clear('list');
end