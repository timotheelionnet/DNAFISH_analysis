function flist = get_clean_file_list_prefixsuffix(...
    inputdir,inc,exc,inc_pref,inc_suf, exc_pref,exc_suf,recursive,case_sensitive)
%input: a directory path (inputdir)
%a cell array containing strings of characters that are included in the
%files to be kept (inc) - set to empty if unnecessary
%a cell array containing strings of characters that are included in the
%files to be excluded (exc)- set to empty if unnecessary
%a number (0 or 1) that sets the recursive option (1 = recursive; 0 = root
%folder only)
%case sensitive: set to 0 or 1
%outputs a list of all the file names (full file name including path)

%% collect all files in folder, including subfolders if recursive option
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
   numel(flist);
end

%% enforcing the file names (including path) includes the <inc> strings
if ~isempty(inc)
    for i=1:numel(inc)
        tmpinc = inc{i};
        if ~isempty(inc{i})
            if ~case_sensitive
                x = regexpi(flist,tmpinc);
            else
                x = regexp(flist,tmpinc);
            end
            idx = ~(cellfun(@isempty,x));
            flist = flist(logical(idx));
        end
    end  
end

%% enforcing the file names (including path) DONT include the <exc> strings
if ~isempty(exc)
    for i=1:numel(exc)
        tmpexc = exc{i};
        if ~isempty(exc{i})
            if ~case_sensitive
                x = regexpi(flist,tmpexc);
            else
                x = regexp(flist,tmpexc);
            end
            idx = cellfun(@isempty,x);
            flist = flist(logical(idx));
        end
    end  
end
if isrow(flist)
    flist = flist';
end

%% enforcing the file names have the right prefix
[~,flist2] = cellfun(@fileparts,flist,'UniformOutput',0);

if ~isempty(inc_pref)
    for i=1:numel(inc_pref)
        tmpinc = inc_pref{i};
        if ~isempty(inc_pref{i})
            if ~case_sensitive
                x = regexpi(flist2,tmpinc);
            else
                x = regexp(flist2,tmpinc);
            end
            
            y = zeros(size(x));
            for j=1:numel(x)
                if isempty(x{j})
                    y(j) = 0;
                else
                    y(j) = min(x{j});
                end
            end
            flist = flist(y==1);
        end
    end  
end

%% enforcing the file names have the right suffix
%redoing this at each step because flist might have been ammended

[~,flist2] = cellfun(@fileparts,flist,'UniformOutput',0);

if ~isempty(inc_suf)
    for i=1:numel(inc_suf)
        tmpinc = inc_suf{i};
        if ~isempty(inc_suf{i})
            if ~case_sensitive
                [~,x] = regexpi(flist2,tmpinc);
            else
                [~,x] = regexp(flist2,tmpinc);
            end
            
            y = zeros(size(x));
            for j=1:numel(x)
                if isempty(x{j})
                    y(j) = 1;
                else
                    y(j) = max(x{j}) - numel(flist2{j});
                end
            end
            flist = flist(y==0);
        end
    end  
end

%% enforcing the file names DONT have the excluded prefix
[~,flist2] = cellfun(@fileparts,flist,'UniformOutput',0);

if ~isempty(exc_pref)
    for i=1:numel(exc_pref)
        tmpexc = exc_pref{i};
        if ~isempty(exc_pref{i})
            if ~case_sensitive
                x = regexpi(flist2,tmpexc);
            else
                x = regexp(flist2,tmpexc);
            end
            
            y = zeros(size(x));
            for j=1:numel(x)
                if isempty(x{j})
                    y(j) = 0;
                else
                    y(j) = min(x{j});
                end
            end
            flist = flist(y~=1);
        end
    end  
end
%% enforcing the file names DONT have the excluded suffix
%redoing this at each step because flist might have been ammended

[~,flist2] = cellfun(@fileparts,flist,'UniformOutput',0);

if ~isempty(exc_suf)
    for i=1:numel(exc_suf)
        tmpexc = exc_suf{i};
        if ~isempty(exc_suf{i})
            if ~case_sensitive
                [~,x] = regexpi(flist2,tmpexc);
            else
                [~,x] = regexp(flist2,tmpexc);
            end
            
            y = zeros(size(x));
            for j=1:numel(x)
                if isempty(x{j})
                    y(j) = 1;
                else
                    y(j) = max(x{j}) - numel(flist2{j});
                end
            end
            flist = flist(y~=0);
        end
    end  
end

end