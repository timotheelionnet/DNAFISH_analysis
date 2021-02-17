function Idx = find_matching_spots(dr,thresh)

%dr is a matrix n1 x n2 of distances between spots from channel 1 and
%channel 2

%outputs the indices of matching spots Idx = n1 x 1 
%(zero when no match)


%% separate dataset into connected components (isolated group of spots within threshold from each other)
[n1,n2] = size(dr);

Idx_left = 1:n1;
conn = [];
ncon = 0;
while ~isempty(Idx_left)
    ncon = ncon+1;
    conn{ncon} = {Idx_left(end); []};
    Idx_left = Idx_left(1:end-1);
    
    spot_found = 1;
    while spot_found
        spfound2 = [];
        for i = 1:numel(conn{ncon}{1}) 
            k = find(dr(conn{ncon}{1}(i),:) < thresh);
            spfound2 = [spfound2,k];
        end
        spfound2 = unique(spfound2);
        spfound2 = spfound2(~ismember(spfound2,conn{ncon}{2}));
        if isempty(spfound2)
            spot_found = 0;
        else
            conn{ncon}{2} = [conn{ncon}{2},spfound2];
            spfound1 = [];
            for i = 1:numel(conn{ncon}{2}) 
                k = find(dr(:,conn{ncon}{2}(i)) < thresh)';
                spfound1 = [spfound1,k];
            end
            spfound1 = unique(spfound1);
            spfound1 = spfound1(~ismember(spfound1,conn{ncon}{1}));
            
            conn{ncon}{1} = [conn{ncon}{1},spfound1];
        end
    end
    Idx_left = setdiff(Idx_left,conn{ncon}{1});
end

xxx_all1 = [];
xxx_all2 = [];
size_all1 = [];
size_all2 = [];
for i=1:numel(conn)
    xxx_all1 = [xxx_all1,conn{i}{1}];
    xxx_all2 = [xxx_all2,conn{i}{2}];
    size_all1(i) = numel(conn{i}{1});
    size_all2(i) = numel(conn{i}{2});
end


%% go through each connected components and try to optimize the pairing
%approach is brute force, computation time will increase
%dereasonably quite fast if the connected components include too many
%spots.
Idx = [];
for i=1:numel(conn)

    if isempty(conn{i}{2})
        Idx = [Idx;[conn{i}{1},0,NaN]];
    else
        m1 = numel(conn{i}{1});
        m2 = numel(conn{i}{2});
        m = min(m1,m2);
        
        i1 = repmat(conn{i}{1}',1,m2);
        i2 = repmat(conn{i}{2},m1,1);
        l = sub2ind(size(dr),i1,i2);
        drtmp = dr(l);
        
        p = find(drtmp < thresh);
        p = reshape(p,numel(p),1);
        [p1,p2] = ind2sub(size(drtmp),p);
        drtmp = drtmp(p);
        drtmp = reshape(drtmp,numel(drtmp),1);
        p = [p1,p2,drtmp];
        
        pmax = 0;
        mmax = m;
        max_achieved = 0;
        best_pair = [];
        while mmax > 0 && ~max_achieved
            ppairs = nchoosek(1:size(p2,1),mmax);
            j=1;
            while j<=size(ppairs,1) && ~max_achieved
                ptmp = p(ppairs(j,:),1:3);
                if numel(unique(ptmp(:,1))) == mmax && numel(unique(ptmp(:,2))) == mmax 
                    max_achieved = 1;
                    best_pair = ptmp;
                    jmax = j;
                end
                j=j+1;
            end
            mmax = mmax - 1;
        end
        
        if ~isempty(best_pair)
            best_pair(:,1) = i1(best_pair(:,1),1);
            best_pair(:,2) = i2(1,best_pair(:,2))';
            
            conn{i}{3} = best_pair;
            leftover1 = setdiff(conn{i}{1},best_pair(:,1))';
            if ~isempty(leftover1)
                leftover1 = [leftover1,zeros(size(leftover1,1),1),NaN(size(leftover1,1),1)];
            end
            
            leftover2 = setdiff(conn{i}{2},best_pair(:,2))';
            if ~isempty(leftover2)
                leftover2 = [zeros(size(leftover2,1),1),leftover2,NaN(size(leftover2,1),1)];
            end
            Idx = [Idx;best_pair;leftover1;leftover2];
        end
    end
end

leftover2 = setdiff(1:n2,unique(Idx(:,2)))';
if ~isempty(leftover2)
    leftover2 = [zeros(size(leftover2,1),1),leftover2,NaN(size(leftover2,1),1)];
end
Idx = [Idx;leftover2];
Idx = sortrows(Idx);
Idx = [Idx(Idx(:,1)~=0,:);Idx(Idx(:,1)==0,:)];

end
