function [Idx,dt] = find_acq_script_from_MM_metadata2( metadata_fname)
%Idx: info about the Micromanager tiff slices, sotred into a 3 column array (frame / Channel / Z)
%each row gives the corresponding TCZ values of the corresponding slice in
%the tiff stack.

%dt: average interval between frames in seconds measured from the time stamps

%load metadata
fid=fopen(metadata_fname,'rt');
data=textscan(fid,'%s');
fclose(fid);

data = data{1,1};

%collect the frame # - channel - z info from the metadata
xFrameKey = cellfun(@strfind,data,repmat({'"FrameKey'},size(data)),'UniformOutput',0);

xFrameKey = ~cellfun(@isempty,xFrameKey); 

idata = data(xFrameKey);

idata = cellfun(@strrep,idata,repmat({'"FrameKey-'},size(idata)), repmat({''},size(idata)),'UniformOutput',0);

idata = cellfun(@strrep,idata,repmat({'":'},size(idata)), repmat({''},size(idata)),'UniformOutput',0);

idata = cellfun(@strrep,idata,repmat({'-'},size(idata)), repmat({' '},size(idata)),'UniformOutput',0);

idata = char(idata); 

%store info into a 3 column array (frame / Channel / Z)
for i =1:size(idata,1)
    Idx(i,1:3) = sscanf(idata(i,:),'%f %f %f');
end



%collect the frame times

xStartms = cellfun(@strfind,data,repmat({'"Time":'},size(data)),'UniformOutput',0);

xStartms = ~cellfun(@isempty,xStartms);

tdata = data([false;false;xStartms(1:end-2,1)]);

tdata = cellfun(@strrep,tdata,repmat({':'},size(tdata)), repmat({' '},size(tdata)),'UniformOutput',0);

tdata = char(tdata); 

%store info into a 3 column array (hour / min / s)
for i =1:size(tdata,1)
    t(i,1:3) = sscanf(tdata(i,:),'%f %f %f');
end

t = t(:,1)*3600 + t(:,2)*60 + t(:,3);

timepoints_counted = numel(tdata)
frames_IDed = max(Idx(:,1))+1

dt = (max(t)- min(t))/(frames_IDed-1)

end
