function res = is_locfile_format(fname)
% this functions checks that a file has the expected format

    if ischar(fname) == 0
        res = 0;
        return;
    else
        [~,~,ext] = fileparts(fname);
        
        %%modify the next line to add a recognized format
        res = strcmp(ext,'.loc') || strcmp(ext,'.loc3'); 
    end
    
end