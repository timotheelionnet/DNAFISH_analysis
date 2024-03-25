function maxproj_tif_filelist( data_flist,save_folder )
%given a Micromanager file pair metadata_fname,data_fname (or list, data_flist in the case MM saved one large acuiqisiotn in mutliple tiffs):

%uses the metadata file to reconstruct the correct data stack and maxproject it at each time point, for each color channel.
% 
%generate one movie stack for each color channel
%
% Also saves log files that tell what the conversion did, 
%time index files for each channel so that each frame of the max projected movie can be related to its corresponding acquisition time.

%Note that if stack #1 (say, A56_003_MMStack_Pos0.ome.tif) contains 2673 cycles plus one frame, 
%the program will split it into its color channels for the first 2673 cycles. 
%The remaining frame is added at the beginning of stack #2 (say, A56_003_MMStack_Pos0_1.ome.tif), 
%so each channel stack starts at the beginning of an imaging cycle.

%%
if isempty(save_folder)
    save_in_same_folder = 1;
else
    save_in_same_folder = 0;
    if exist(save_folder,'file') ~= 7
        mkdir(save_folder);
    end
end 

%start data log
NL = sprintf('\r\n');
file_str = ['TIF 3D max-proj Conversion Log',NL,NL];
tmp_str = [];


%for i=2:numel(data_flist)
for i=1:numel(data_flist)
    i
    if ~isempty(data_flist{i})
        %load data; concatenate with remaining unfinished cycles from previous
        %imaging stack if needed.
        data_flist{i}
        data = tiffread6(data_flist{i});
        [dirname,fname,~] = fileparts(data_flist{i});
        if save_in_same_folder
            save_folder = dirname;
        end
        
        %log the data file name
        tmp_str = [tmp_str ,'Data File #',num2str(i),': ',data_flist{i},NL,NL];   

        if isempty(data)
            tmp_str = [tmp_str ,'Could not load data',NL,NL];
        else
            
            %find out data size
            [nx,ny,slice_end] = size(data)
            tmp_str = [tmp_str ,'size',num2str(nx),' ',num2str(ny),' ',num2str(slice_end),NL,NL];   
            
            if nx == 0 || ny == 0 || slice_end == 0
                tmp_str = [tmp_str ,'Could not load data',NL,NL];
            else
                C = squeeze(max(data,[],3));
                size(C)

                %save entire stack
                disp('saving tiff...');
                disp(fullfile(save_folder,[fname,'_max.tif']));
                saveastiff2(C,fullfile(save_folder,[fname,'_max.tif']));

                %disp('saving avi...');
                %save_stack_as_avi(C,fullfile(dirname,[fname,'_C',num2str(j),'.avi']),'slice');
            end
        end
    end
end

file_str = [file_str, tmp_str];

[dirname,fname,~] = fileparts(data_flist{1});    
fid1 = fopen(fullfile(dirname,[fname,'_Conversion_log.txt']),'w');
fprintf(fid1,'%s \r\n',file_str);
fclose(fid1);

clear('data');
end

