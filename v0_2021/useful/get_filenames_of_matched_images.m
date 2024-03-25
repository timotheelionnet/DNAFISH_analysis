function fnames = get_filenames_of_matched_images()
%filenames should have the following structure:
%base_name_channel1string_randomstring.tif
%basename should be the same bewteen the three color channels    
    
persistent c1_dir
if isempty(c1_dir) 
    c1_dir = '\'; 
end


%% user selects filenames and directories
[c1_fname,c1_dir,fidx] = uigetfile('*.tif','Chose First Red Channel Image',c1_dir);
if fidx == 0
    return
end

[c2_fname,c2_dir,fidx] = uigetfile('*.tif','Chose First Green Channel Image',c1_dir);
if fidx == 0
    return
end

[c3_fname,c3_dir,fidx] = uigetfile('*.tif','(Optional) Chose First Blue Channel Image',c1_dir);
if fidx == 0
    numchannels = 2;
else
    numchannels = 3;
end

%% determine base name and extensions
% extract the base name 
if numchannels == 2
    %dummy 3rd channel to simply following code
    c3_fname = c2_fname;
end

[~,~,ext1,ext2,ext3] = find_basename_and_extensions(c1_fname,c2_fname,c3_fname);



%% establish the list of matched filenames

list1 = dir(c1_dir);
list2 = dir(c2_dir);
if numchannels == 3
    list3 = dir(c3_dir);
end


fnames = cell(1,numchannels);
nimgs = 0;

for i=1:length(list1) 
    if ~ list1(i).isdir && is_recognized_img_format(list1(i).name)%checking that it is an image file
        if ~isempty( strfind(list1(i).name, ext1 ))
            nimgs = nimgs+1;
            fnames{nimgs,1} = [c1_dir,list1(i).name];
            [prefix,suffix] = find_basename(list1(i).name,ext1);
            
            j=1;
            while j<=length(list2)
                if ~isempty( strfind(list2(j).name, prefix )) && ~isempty( strfind(list2(j).name, suffix )) && ~isempty( strfind(list2(j).name, ext2 ))
                    fnames{nimgs,2} = [c2_dir,list2(j).name];
                    break
                else
                    j = j+1;
                end
            end
            
            if numchannels == 3
                k=1;
                while k<=length(list3)
                    if ~isempty( strfind(list3(k).name, prefix )) && ~isempty( strfind(list3(k).name, suffix )) && ~isempty( strfind(list3(k).name, ext3 ))
                        fnames{nimgs,3} = [c3_dir,list3(k).name];
                        break
                    else
                        k = k+1;
                    end
                end
            end

        end
        
    end
    
end



end

function [prefix,suffix,ext1,ext2,ext3] = find_basename_and_extensions(c1_fname,c2_fname,c3_fname)

%identifies the common strings at the beginning and end of three filenames:
%example c001_red_xxx.tif, c001_green_xxx.tif, c001_blue_xxx.tif will
%return
%prefix = 'c001_'
%suffix = '_xxx.tif'


%min size of the full name
nchar = min([length(c1_fname),length(c2_fname),length(c3_fname)]);

%find a potential common string at the start of the filename
i=1;
while i <=nchar
    res = strncmp(c1_fname,c2_fname,i);
    res = res*strncmp(c1_fname,c3_fname,i);
    if res == 0 
        break 
    else
        i = i+1;
    end    
end

if i == 1
    prefix = '';
elseif i == nchar
    disp('could not identify channel specific string');
    return
else    
    prefix = c1_fname(1:i-1);
end

%find a potential common string at the end of the filename
ic1_fname = flipud(c1_fname')';
ic2_fname = flipud(c2_fname')';
ic3_fname = flipud(c3_fname')';

j=1;
while j <=nchar
    res = strncmp(ic1_fname,ic2_fname,j);
    res = res*strncmp(ic1_fname,ic3_fname,j);
    if res == 0 
        break 
    else
        j = j+1;
    end    
end

if j==1
    suffix = '';    
else
    suffix = flipud(ic1_fname(1:j-1)')';
end

if i == 1 && j ==1
    disp('could not identify basename');
    ext1 = '';
    ext2 = '';
    ext3 = '';
    return
end

% extract the channel-specific extensions
l1 = length(prefix);
l2 = length(suffix);


ext1 = c1_fname(l1+1:length(c1_fname)-l2);
ext2 = c2_fname(l1+1:length(c2_fname)-l2);
ext3 = c3_fname(l1+1:length(c3_fname)-l2);

end

function [prefix,suffix] = find_basename(fname,ext)

    k = strfind(fname,ext);
    prefix = fname(1:k-1);
    suffix = fname(k+length(ext):length(fname));

end


