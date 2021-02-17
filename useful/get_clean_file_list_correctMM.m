function flist = get_clean_file_list_correctMM(inputdir,inc, exc,recursive,case_sensitive)
%input: a directory path (inputdir)
%a cell array containing strings of characters that are included in the
%files to be kept (inc) - set to empty if unnecessary
%a cell array containing strings of characters that are included in the
%files to be excluded (exc)- set to empty if unnecessary
%a number (0 or 1) that sets the recursive option (1 = recursive; 0 = root
%folder only)
%case sensitive: set to 0 or 1
%outputs a list of all the file names (full file name including path)

if recursive == 0
   flist = dir( inputdir);
   %remove directories
   flist = flist(~[flist.isdir]);
   dircell = cell(size(flist));
   dircell(:) = {inputdir};
   fnames = struct2cell(flist);
   fnames = fnames(1,:);
   fnames = fnames';
   flist = cellfun(@fullfile,dircell,fnames,'UniformOutput',0);
else
    %get all files
   flist = dirrec(inputdir);
end

%create mirror list where the MM wavelengths strings are placed at the end
flist2 = flist;
for i=1:numel(flist)
    k1 = strfind(flist{i},'_w');
    k2 = strfind(flist{i},'Cube');
    k3 = strfind(flist{i},'_');
    k4 = strfind(flist{i},'.');
    if ~isempty(k4)
        k4 = k4(end);
    else
        k4 = length(flist{i});
    end
    if ~isempty(k1) && ~isempty(k2) && ~isempty(k3)
        k3 = k3(k3>k2);
        k3 = min(k3);
        flist{i} = [flist{i}(1:k1-1),flist{i}(k3:k4-1),flist{i}(k1:k3-1),flist{i}(k4:end)];
    end
end

if ~isempty(inc)
    for i=1:numel(inc)
        tmpinc = inc{i};
        if ~isempty(inc{i})
            if ~case_sensitive
                x = regexpi(flist,tmpinc);
            else
                x = regexp(flist,tmpinc);
            end
            idxinc = ~(cellfun(@isempty,x));
        else
            idxinc = ones(size(flist));
        end
    end
else
    idxinc = ones(size(flist));
end


if ~isempty(exc)
    for i=1:numel(exc)
        tmpexc = exc{i};
        if ~isempty(exc{i})
            if ~case_sensitive
                x = regexpi(flist,tmpexc);
            else
                x = regexp(flist,tmpexc);
            end
            idxexc = cellfun(@isempty,x);
        else
            idxexc = ones(size(flist));
        end
    end 
else
    idxexc = ones(size(flist));
end

%combine inclusion and exclusion condition
idx = logical(idxinc.*idxexc);

%retrieve original filenames that satisfy the condition
flist = flist2(idx);

if isrow(flist)
    flist = flist';
end

end