function [fname,basename] = assemble_multiple_file_lists_correctMM(listing,suffix)

listing2 = listing;

%move wavelength string in MM name to the end
for i=1:numel(listing)
    for j=1:size(listing{i},1)
        k1 = strfind(listing{i}{j},'_w');
        k2 = strfind(listing{i}{j},'Cube');
        k3 = strfind(listing{i}{j},'_');
        k4 = strfind(listing{i}{j},'.');
        if ~isempty(k4)
            k4 = k4(end);
        else
            k4 = length(flist{i});
        end
        if ~isempty(k1) && ~isempty(k2) && ~isempty(k3)
            k3 = k3(k3>k2);
            k3 = min(k3);
            listing{i}{j}= [listing{i}{j}(1:k1-1),listing{i}{j}(k3:k4-1),listing{i}{j}(k1:k3-1),listing{i}{j}(k4:end)];
        end
    end    
end

numfiles = 0;
basename = cell('');
for i1=1:size(listing{1},1)
    %search for candidates in the channel 1
    [~,basenametmp,~] = fileparts(listing{1}{i1});
    k1 = strfind(basenametmp,suffix{1});
    basenametmp = basenametmp(1:k1-1);
    idx(1) = i1;

    %once candidate file is found in channel 1,
    %look into other directories for files with matching basename and the
    %adequate suffix
    find_all_files =1;
    for i = 2:numel(listing)
        flag = 0;
        for i2 = 1:size(listing{i},1)
            [~,name2,~] = fileparts(listing{i}{i2});
            k2 = strfind(name2,basenametmp);
            if isscalar(k2) 
                flag = 1;
                idx(i) = i2;
                break;
            end
        end 
        find_all_files = find_all_files*flag;
    end

    %store file names if all matching names have been found
    if find_all_files == 1
        numfiles = numfiles +1;
        for i = 1:numel(listing)
            fname{numfiles,i} = fullfile(listing2{i}{idx(i)});
            basename{numfiles} = basenametmp;
        end
    end

end

end






