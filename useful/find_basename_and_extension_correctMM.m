function [path,basename, prefix,suffix,ext] = find_basename_and_extension_correctMM(varargin)

%load and parse file names
for i=1:nargin
    [path{i},name{i},ext{i}] = fileparts(varargin{i}); 
end

%move MetaMorph wavelength string towards the end of each filename
for i=1:numel(name)
    k1 = strfind(name{i},'_w');
    k2 = strfind(name{i},'Cube');
    k3 = strfind(name{i},'_');
    if ~isempty(k1) && ~isempty(k2) && ~isempty(k3)
        k3 = k3(k3>k2);
        k3 = min(k3);
        name{i} = [name{i}(1:k1-1),name{i}(k3:end),name{i}(k1:k3-1)];
    end
    l(i) = length(name{i});
end

%find shared elements at the beginning and the end of the file
lmin = min(l);%smallest name length
i1 = 0; %if there is a shared substring at the beginining of the file, i1 is the last index of the substring
i2 = 0; %if there is a shared substring at the end of the file, i2 is the index of the first character of that substring
for i=1:lmin
    is_string_start_similar = 1;
    for j = 2:nargin
        is_string_start_similar = is_string_start_similar * strcmp(name{1}(1:i),name{j}(1:i)) ;  
    end
    
    if is_string_start_similar == 1
        i1 = i;
    end
    
    is_string_end_similar = 1;
    for j = 2:nargin
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
    for j = 1:nargin
        suffix{j} = name{j}(i1+1:l(j));
    end
elseif i1 == 0 && i2 > 0
    %common substring after a variable prefix
    basename = name{1}(l(1)-i2+1:l(1));
    for j = 1:nargin
        prefix = name{j}(1:l(j)-i2);
    end
    suffix = {};
    
elseif i1 > 0 && i2 > 0
    %beginning and end are identical; variable part (stored in suffix) comes in between
    %basename is now a structure with 2 entries
    basename = {name{1}(1:i1),name{1}(l(1)-i2+1:l(1))};
    prefix = {};
    for j = 1:nargin
        suffix{j} = name{j}(i1+1:l(j)-i2);
    end
end
