
addpath("Dependencies\CensoredRepeatSequences\");

human_repeats = fastaread('./Dependencies/CensoredRepeatSequences/HumanRepeats.fasta');
%% FC FISH Probe design
% 7/10/2023

% Input fasta for region of interest

filePath = "G:\Finn\20230710_DEVIN_ECDNA_FISH\hg38_Devin_enhancer.fa";
region = fastaread(filePath);
[saveDir,~] = fileparts(filePath);

% split into probes of length k
k = 33;

% gc range
gc_range = [0.2, 0.8];

% tm min
my_Tm_min =42;

% look for homology to repeats
repeat_homol = 14;

% filter for intra set homology
set_homol = 12;

%%
% get seq
seq = region.Sequence;
header = region.Header;


seq_len = length(seq);

indx = 1:k:seq_len;

% loop thru the seq and get kmers
all_kmers = [];
for i = 1:numel(indx)
    
    % check if we have less than k bases left
    if indx(i) + k-1> seq_len
        continue
    end

    % get the kmer
    my_kmer = seq(indx(i):indx(i)+ k-1);

    % add to list
    all_kmers = vertcat(all_kmers, my_kmer);

end

disp('Starting with ' + string(height(all_kmers)) + ' probe candidates');

%% load repeats

% human_repeats = fastaread("G:\Finn\20230710_DEVIN_ECDNA_FISH\CensoredRepeatSequences\HumanRepeats.fasta");

% combine them into one mega string!!!!

forbidden_seq = upper([human_repeats.Sequence]);




%% start filtering kmers
gc_rmv_cnt = 0;
tm_rmv_cnt = 0;
repeat_rmv_cnt = 0;

% store gc for plotting
my_gc_dist =[];


my_Tm_dist = [];

% search thru unfiltered kmers
filtered_kmers = [];

for i = 1:height(all_kmers)
    
    % cur kmer
    my_kmer = all_kmers(i, :);

    % get gc and filter ---------------------------
    cur_gc = get_gc_content(my_kmer);

    if cur_gc < gc_range(1) || cur_gc > gc_range(2)

        gc_rmv_cnt = gc_rmv_cnt +1;
        continue
    else
        my_gc_dist = [my_gc_dist cur_gc];

    end

    % check melt temp and filter------------------

    cur_kmer_prop = oligoprop(my_kmer);

%     cur_kmer_prop.Tm

    cur_Tm = max(cur_kmer_prop.Tm); 
    my_Tm_dist = [my_Tm_dist cur_Tm];

    if cur_Tm < my_Tm_min
        
        tm_rmv_cnt = tm_rmv_cnt +1;
        continue
    else
        

    end

    
    kmer_len = length(my_kmer); 
    % break into 14 nt submers
    submer_indx = 1:repeat_homol-1:kmer_len;
    
    % reset for each kmer
    forbidden_kmer = 0;
    % loop thru submers and look for homology
    for j = 1:numel(submer_indx)

        if forbidden_kmer == 1
            continue
        end
    
        % check if we have less than k bases left
        if submer_indx(j) + repeat_homol-1 > kmer_len
            continue
        end
    
        % get the submer
        my_submer = my_kmer(submer_indx(j):submer_indx(j)+ repeat_homol-1);
        
        % does it have homology to forbidden repeat seqs?
        is_forbidden = strfind(forbidden_seq, my_submer);

        if ~isempty(is_forbidden)
            forbidden_kmer = 1;
            
        end


        
        
    

    end

    if forbidden_kmer == 1
        repeat_rmv_cnt = repeat_rmv_cnt +1;
        
        continue
    end

    % add surviving kmer to list of filtered kmers
    filtered_kmers = vertcat(filtered_kmers, my_kmer);


end

%%

filtered_kmers2 = [];
self_homol_cnt = 0;

for i = 1:height(filtered_kmers)
    % get ith kmer
    my_kmer = filtered_kmers(i, :);

    kmer_len = length(my_kmer); 

    % break into specificed len submers
    submer_indx = 1:set_homol-1:kmer_len;
    
    % switch for if kmer is disqualified 
    forbidden_kmer = 0;

    % loop thru submers and look for homology to all other kmers
    % cnt how intra set hits for each submer
    set_homol_cnt = 0;
    for j = 1:numel(submer_indx)


        
        % dont waste time, if failed once, skip
        if forbidden_kmer == 1
            continue
        end
    
        % check if we have less than set_homol bases left
        if submer_indx(j) + set_homol-1 > kmer_len
            continue
        end
    
        % get the submer
        my_submer = my_kmer(submer_indx(j):submer_indx(j)+ set_homol-1);
        
        % does it have homology to forbidden repeat seqs?
        % for current submer j, loop thru all k kmers and look for homology
        for k = 1:height(filtered_kmers)

            if k == i % skip self
                continue
            end

            query_kmer = filtered_kmers(k,:);

            % does current submer have homology to another probe
            is_forbidden = strfind(query_kmer, my_submer); % does submer align?
            

            if ~isempty(is_forbidden)
                % count +1 for each probe the submer has homology to
                set_homol_cnt = set_homol_cnt +1;
                
                
            end
        end

        if set_homol_cnt >= 1 % homologous to self and one other
            forbidden_kmer = 1;
        end

    

    end
    
    % count how many probes have been removed due to having excess homology
    % to other probes in the set
    if forbidden_kmer == 1
        self_homol_cnt = self_homol_cnt +1;
        
        continue
    end

    % add surviving kmer to list of filtered kmers
    filtered_kmers2 = vertcat(filtered_kmers2, my_kmer);
end


%%

% -------------------------------------------------------
% report filtering results
% ------------------------------------------------------
disp('~~~~~~~~~~~~~~~~~~~~');
disp('Started with ' + string(height(all_kmers)) + ' probe candidates');
disp('gc range: ' + string(gc_range(1)) + '-' + string(gc_range(2)) );
disp('probes removed due to gc: ' + string(gc_rmv_cnt))

disp('Tm thresh (C): ' + string(my_Tm_min))
disp('probes removed due to Tm: ' + string(tm_rmv_cnt))

disp('repeat homology thresh: ' + string(repeat_homol))
disp('probes removed due to repeat homology: ' + string(repeat_rmv_cnt))

disp('inter probe homology thresh: ' + string(set_homol))
disp('probes removed due to inter-probe homology: ' + string(set_homol_cnt))

disp('Remaining Probes: ' + string(height(filtered_kmers)))

disp('~~~~~~~~~~~~~~~~~~~~');

% save params to txt
param_str = ['Started with ' + string(height(all_kmers)) + ' probe candidates'...,
            'gc range: ' + string(gc_range(1)) + '-' + string(gc_range(2)) ...,
            'probes removed due to gc: ' + string(gc_rmv_cnt)...,
            'Tm thresh (C): ' + string(my_Tm_min)...,
            'probes removed due to Tm: ' + string(tm_rmv_cnt)...,
            'repeat homology thresh: ' + string(repeat_homol)...,
            'probes removed due to repeat homology: ' + string(repeat_rmv_cnt)...,
            'inter probe homology thresh: ' + string(set_homol)...,
            'probes removed due to inter-probe homology: ' + string(set_homol_cnt)...,
            'Remaining Probes: ' + string(height(filtered_kmers))            
            ]';


%% save results

savename = fullfile(saveDir, 'probes.csv');
writematrix(filtered_kmers2, savename);

params_savename = fullfile(saveDir, 'design_params.txt');
writematrix(param_str, params_savename);


disp('saved to')
disp(savename)

%% instructions from boettiger and zhuang labs

% For each block, pick non-overlapping 40-nt sequences. 
% These 40-nt sequences must have (i) 
% a melting temperature >65 °C, (ii) 
% a G/C range of 20–80%, (iii) 
% a homology <12 nt to any other barcode in the pool and (iv) 
% a homology <14 nt to repetitive elements in the genome.