
function [sp1,sp2] = sort_loc_data_nearest_neighbors(r1,r2)

%% first sort nearest neighbor pairs of spots from each input loc file
%store in sp1, and sp2 formated as
%[x1,x2,x3,...;
% y1,y2,y3,...]

frames = unique(r1(:,end));
nframes = numel(frames);

sp1 = [];
sp2 = [];
comp1 = [];
comp2 = [];

for i=1:nframes
    r1tmp = r1(r1(:,end) == frames(i),:);
    r2tmp = r2(r2(:,end) == frames(i),:);   
    
    [r1tmp,r2tmp] = assign_spot_pairs(r1tmp,r2tmp,0,0); 
    
    sp1 = [sp1, [  r1tmp(:,1:2)' ; ones(1,size(r1tmp,1)) ]  ];
    sp2 = [sp2, [  r2tmp(:,1:2)' ; ones(1,size(r2tmp,1)) ]  ];
    
    %save the remaining columns (e.g. Intensity, frame) so we can add them back to the corrected
    %spot position columns at the end
    if size(r1,2) >2 && size(r2,2) >2 
        comp1 = [comp1; r1tmp(:,3:end)];
        comp2 = [comp2; r2tmp(:,3:end)];
    end
    
end

sp1 = [sp1', comp1];
sp2 = [sp2', comp2];
end
