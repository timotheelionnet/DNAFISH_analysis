function Idx = find_acq_script_from_MM_metadata( metadata_fname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%collect the time - channel - z info from the metadata
fid=fopen(metadata_fname,'rt');
data=textscan(fid,'%s');
fclose(fid);

data = data{1,1};

xFrameKey = cellfun(@strfind,data,repmat({'"FrameKey'},size(data)),'UniformOutput',0);

xFrameKey = ~cellfun(@isempty,xFrameKey); 

data = data(xFrameKey);

data = cellfun(@strrep,data,repmat({'"FrameKey-'},size(data)), repmat({''},size(data)),'UniformOutput',0);

data = cellfun(@strrep,data,repmat({'":'},size(data)), repmat({''},size(data)),'UniformOutput',0);

data = cellfun(@strrep,data,repmat({'-'},size(data)), repmat({' '},size(data)),'UniformOutput',0);

data = char(data); 

for i =1:size(data,1)
    Idx(i,1:3) = sscanf(data(i,:),'%f %f %f');
end

%open tif file

end

