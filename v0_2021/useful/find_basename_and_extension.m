function [path,basename, prefix,suffix,ext] = find_basename_and_extension(varargin)

%load and parse file names
k=1;
for i=1:nargin
    if ~iscell(varargin{i})
        [path{k},name{k},ext{k}] = fileparts(varargin{i}); 
        l(k) = length(name{k});
        k = k+1;
    else
        for j=1:numel(varargin{i})
            [path{k},name{k},ext{k}] = fileparts(varargin{i}{j}); 
            l(k) = length(name{k});
            k = k+1;
        end
    end
end
nfiles = k-1;

%find shared elements at the beginning and the end of the file
lmin = min(l);%smallest name length
i1 = 0; %if there is a shared substring at the beginining of the file, i1 is the last index of the substring
i2 = 0; %if there is a shared substring at the end of the file, i2 is the index of the first character of that substring
is_string_start_similar = 1;
is_string_end_similar = 1;
for i=1:lmin
    for j = 2:nfiles
        is_string_start_similar = is_string_start_similar * strcmp(name{1}(1:i),name{j}(1:i)) ;  
    end
    
    if is_string_start_similar == 1
        i1 = i;
    end
    
    for j = 2:nfiles
        is_string_end_similar = is_string_end_similar * strcmp(name{1}(l(1)-i+1:l(1)),name{j}(l(j)-i+1:l(j))) ;  
    end
    
    if is_string_end_similar
        i2 = i;
    end
end


%going through the different cases
if i1 == 0 && i2 == 0
    %no shared substring found at beginning or end
    basename = [];
    suffix = name;
    prefix = {};
    
elseif i1 > 0 && i2 == 0
    %common substring followed by a variable suffix
    basename = name{1}(1:i1);
    prefix = {};
    for j = 1:nfiles
        suffix{j} = name{j}(i1+1:l(j));
    end
elseif i1 == 0 && i2 > 0
    %common substring after a variable prefix
    basename = name{1}(l(1)-i2+1:l(1));
    for j = 1:nfiles
        prefix = name{j}(1:l(j)-i2);
    end
    suffix = {};
    
elseif i1 > 0 && i2 > 0
    %beginning and end are identical; variable part (stored in suffix) comes in between
    %basename is now a structure with 2 entries
    basename = {name{1}(1:i1),name{1}(l(1)-i2+1:l(1))};
    prefix = {};
    for j = 1:nfiles
        suffix{j} = name{j}(i1+1:l(j)-i2);
    end
end
