function movieInfo = output_trajectories_as_structure_from_res2(res)


%formatting for the u-track program suite
%For a movie with N frames, movieInfo is a structure array with N entries.
%Every entry has the fields xCoord, yCoord, zCoord (if 3D) and amp.
%If there are M features in frame i, each one of these fields in
%moveiInfo(i) will be an Mx2 array, where the first column is the value
%(e.g. x-coordinate in xCoord and amplitude in amp) and the second column
%is the standard deviation. If the uncertainty is unknown, make the second
%column all zero.

%the imput has the format [x, y, z, Int, Im. number]

movieInfo=struct('xCoord',{},'yCoord',{},'zCoord',{},'amp',{});

nframes = numel(res);
for i=1:nframes
    
    npts = size(res{i}.final_pix,1);
    movieInfo(i).xCoord = zeros(npts,2);
    movieInfo(i).yCoord = zeros(npts,2);
    movieInfo(i).zCoord = zeros(npts,2);
    movieInfo(i).amp = zeros(npts,2);
    
    movieInfo(i).xCoord(1:npts,1) = res{i}.final_pix(1:npts,1);
    movieInfo(i).yCoord(1:npts,1) = res{i}.final_pix(1:npts,2);
    movieInfo(i).zCoord(1:npts,1) = res{i}.final_pix(1:npts,3);
    movieInfo(i).amp(1:npts,1) = res{i}.final_pix(1:npts,4);
end


end