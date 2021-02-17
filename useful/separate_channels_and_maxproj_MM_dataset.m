function separate_channels_and_maxproj_MM_dataset( metadata_fname,data_fname )

%Find the time - channel - z info relative to each slice from the metadata
Idx = find_acq_script_from_MM_metadata( metadata_fname);

%load data stack
data = tiffread5(data_fname);

%separate channels and save separately
Cidx = unique(Idx(:,2));
[dirname,fname,~] = fileparts(data_fname);

NL = sprintf('\r\n');

file_str = ['Micro-Manager 4D to TIF 3D (z max-projected) Conversion Log',NL,NL,...
    'Data File: ',data_fname,NL,NL,...
    'MetaData File: ',metadata_fname,NL,NL];
    
for i=1:numel(Cidx)
    
    C = data(:,:,Idx(:,2)==Cidx(i));
    [nx,ny,n3] = size(C);
    zslices = numel(unique(Idx(Idx(:,2)==Cidx(i),3)));
    numtimepoints = numel(unique(Idx(Idx(:,2)==Cidx(i),1)));
    
    timepoints = Idx(Idx(:,2)==Cidx(i),:);
    timepoints = timepoints(timepoints(:,3)==max(timepoints(:,3)),1);
    tmp_str = ['Channel ',num2str(i),NL,...
        num2str(n3),' total frames',NL,...
        num2str(zslices),' z-slices',NL];
    
    %if last stack was interrupted
    if zslices*numtimepoints ~= n3
        numtimepoints = numtimepoints - 1;
        C = C(:,:,1:numtimepoints*zslices);
        tmp_str = [tmp_str, num2str(numtimepoints),' time points',NL,...
            num2str(n3 - numtimepoints*zslices),' last frames deleted',NL,NL];
    else
        tmp_str = [tmp_str, num2str(numtimepoints),' time points',NL,NL];
    end
    C = reshape(C,nx,ny,zslices,numtimepoints);
    C = squeeze(max(C,[],3));
    save_as_tiff(C,fullfile(dirname,[fname,'_C',num2str(i),'.tif']));
    clear('C');
    file_str = [file_str, tmp_str];
    save(fullfile(dirname,[fname,'_t',num2str(i),'.txt']),'timepoints','-ascii');
end

fid1 = fopen(fullfile(dirname,[fname,'_Conversion_log.txt']),'w');
fprintf(fid1,'%s \r\n',file_str);
fclose(fid1);
    
save(fullfile(dirname,[fname,'_TCZraw.txt']),'Idx','-ascii');

end

