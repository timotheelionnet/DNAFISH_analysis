function beadchannels_flist = separate_MMchannels_fast_for_beads( metadata_fname,data_flist)
%given a Micromanager file pair metadata_fname,data_fname (or list, data_flist in the case MM saved one large acquisition in mutliple tiffs):

%uses the metadata file to reconstruct the correct data stack (time, color, z) and maxproject it at each time point, for each color channel.
% 
%generate one movie stack for each color channel
%
% Also saves log files that tell what the conversion did, 
%time index files for each channel so that each frame of the max projected movie can be related to its corresponding acquisition time.

%Note that if stack #1 (say, A56_003_MMStack_Pos0.ome.tif) contains 2673 cycles plus one frame, 
%the program will split it into its color channels for the first 2673 cycles. 
%The remaining frame is added at the beginning of stack #2 (say, A56_003_MMStack_Pos0_1.ome.tif), 
%so each channel stack starts at the beginning of an imaging cycle.

%this function is similar to
%separate_channels_and_maxproj_MM_dataset_cleanupframes5,
%but saves less files (e.g. no avi or control tiffs)

%INPUT
%file pair metadata_fname,data_fname 
    %data_fname can be a list (cell array), in the case MM saved one large acquisition in mutliple tiffs

%OUTPUT
%beadchannels_flist{i,j} lists the file name of dataset i, channel j 



%%



%Find the time - channel - z info relative to each slice from the metadata
Idx = find_acq_script_from_MM_metadata( metadata_fname);

%start data log
NL = sprintf('\r\n');
file_str = ['Micro-Manager 4D to TIF 3D (z max-projected) Conversion Log',NL,NL];
file_str = [file_str ,'MetaData File: ',metadata_fname,NL,NL];
tmp_str = [];

%go through all data stacks
slice_start = 0;
data_rem = [];

for i=1:numel(data_flist)
    i
    if ~isempty(data_flist{i})
        %load data; concatenate with remaining unfinished cycles from previous
        %imaging stack if needed.
        data_flist{i}
        size(data_rem)
        data = cat(3,data_rem,tiffread5(data_flist{i}));
        [dirname,fname,~] = fileparts(data_flist{i});

        %log the data file name
        tmp_str = [tmp_str ,'Data File #',num2str(i),': ',data_flist{i},NL,NL];   

        %find out data size
        [nx,ny,slice_end] = size(data);

        %select the indices from the global metadata file that correspond to current data stack
        Idxtmp = Idx(slice_start+1:slice_start+slice_end,:);

        %select an integer number of imaging rounds (including all z and C values). 
        %save the remainder in case the dataset has been truncated in the
        %middle of a z-stack or color acquisition
        Cidx = unique(Idxtmp(:,2));
        nC = numel(Cidx);
        Zidx = unique(Idxtmp(:,3));
        nZ = numel(Zidx);

        %this is the number of full imaging cycles in the dataset (Z + color)
        Ncycles = floor(slice_end/(nC*nZ));
        Idxrem = mod(slice_end,nC*nZ);

        %saving the last unfinished cycle (will be concatenated on top of the
        %following stack)
        if Idxrem>0
            data_rem = data(:,:,end-Idxrem+1:end);
        else
            data_rem = [];
        end

        %setting the start index for the next z-stack (making sure it matches
        %the data with the concatenated unfinished cycle)
        slice_start = slice_start+slice_end - Idxrem;

        %keeping the rest of the data (Ncycles rounds of imaging)
        data = data(:,:,1:end-Idxrem);

        %correcting the current indices
        Idxtmp = Idxtmp(1:Ncycles*nC*nZ,:);

        tmp_str = [tmp_str, num2str(Ncycles),' full Imaging cycles;',NL,...
            num2str(Idxrem),' extra frames from next imaging cycles combined with data from next file.', NL, NL];

        for j=1:numel(Cidx)

            C = data(:,:,Idxtmp(:,2)==Cidx(j));
            [nx,ny,n3] = size(C);
            zslices = numel(unique(Idxtmp(Idxtmp(:,2)==Cidx(j),3)));
            numtimepoints = numel(unique(Idxtmp(Idxtmp(:,2)==Cidx(j),1)));

            timepoints = Idxtmp(Idxtmp(:,2)==Cidx(j),:);
            timepoints = timepoints(timepoints(:,3)==max(timepoints(:,3)),1);
            tmp_str = [tmp_str,'Channel ',num2str(j),NL,...
                num2str(n3),' total frames',NL,...
                num2str(zslices),' z-slices',NL];

            %if last stack was interrupted (overkill sanity check, remnant from earlier version)
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
            size(C)

            %save entire stack
            beadchannels_flist{i,j} = fullfile(dirname,[fname,'_C',num2str(j),'.tif']);
            save_as_tiff(C,beadchannels_flist{i,j});

            save(fullfile(dirname,[fname,'_C',num2str(j),'_time_index.txt']),'timepoints','-ascii');  
            
            clear('C');

        end
    end
end

file_str = [file_str, tmp_str];

[dirname,fname,~] = fileparts(data_flist{1});    
fid1 = fopen(fullfile(dirname,[fname,'_Conversion_log.txt']),'w');
fprintf(fid1,'%s \r\n',file_str);
fclose(fid1);

save(fullfile(dirname,[fname,'_TCZraw.txt']),'Idx','-ascii');

size(data)
clear('data');
end






