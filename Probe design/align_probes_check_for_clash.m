% load pparg locus

% get locus seq
filePath = "G:\Finn\20231122_PPARG_probe_expansion\pparg_prom.fasta";
[saveDir,~] = fileparts(filePath);

locus_seq = fastaread(filePath);
locus_seq = locus_seq.Sequence;

% get probe seqs set 1
probe_seqsT = readtable("G:\Finn\20231122_PPARG_probe_expansion\pparg_promoter_probeset.xlsx");

probe_seqs = probe_seqsT.Var5; 

% get probe seqs set 2
probe_seqsT2 = readtable("G:\Finn\20231122_PPARG_probe_expansion\FC_probes.csv", 'ReadVariableNames', false);

probe_seqs2 = probe_seqsT2.Var1;


%%

probe_aligned = zeros(1, length(locus_seq)); % 0 for every base in the locus seq

probe_indcs = []; % array where we store index of each probe aligned to target

for i = 1:numel(probe_seqs)

    cur_probe = probe_seqs{i};
    
    % find where each oligo lands 
    k = strfind(locus_seq, cur_probe);
    
    probe_indcs = [probe_indcs, k];
end



% now we use probe indcs as indexes in the alignment binary 

probe_aligned(probe_indcs) = 1;

%% visualize where current probe set is bound
figure
plot(1:length(probe_aligned), probe_aligned)

xlabel('location on target (bp)');
title('Set 1')

%% now see where we can fill in the gaps


probe_indcs2 = [];

probe_aligned2 = zeros(1, length(locus_seq)); % 0 for every base in the locus seq

for i = 1:height(probe_seqs2)

    cur_probe = probe_seqs2{i};
    
    % find where each oligo land, %strfind returns starting index of substr
    % within string
    k = strfind(locus_seq, cur_probe);
    
    probe_indcs2 = [probe_indcs2, k];
end

probe_aligned2(probe_indcs2) = 1;

%% visualize where current probe set is bound
figure
plot(1:length(probe_aligned2), probe_aligned2)

%% Figure out where clashing occurs

% each probe in probe set2 must be d bases away

d = 60;

probes2rmv_idx = zeros(size(probe_indcs2));
for i = 1:numel(probe_indcs2)

    cur_idx = probe_indcs2(i);

    % is our probe within distance d of any other probe

    dif_arr = probe_indcs - cur_idx;

    clash_points = abs(dif_arr) <= d;

    if sum(clash_points) >= 1 % does a probe clash 

        probes2rmv_idx(i) = 1; % we will remove the ith probe from the new set
        
    end

end

%% remove clashing probes from probe set 2

disp('~~~~~~~~~~~~~')
disp( sprintf('started with %f probe candidates', numel((probe_seqs2)) ));
disp(strcat(string(sum(probes2rmv_idx)), ' probes in set2 were removed due to clashing'));

disp(sprintf('There remain: %f probes in the new probe set', floor(numel( probe_seqs2(logical(~probes2rmv_idx),1))) ) );

disp('~~~~~~~~~~~~~')

probe_seqs2_filtered = probe_seqs2(logical(~probes2rmv_idx),1);




% align and check 

probe_aligned3 = zeros(1, length(locus_seq)); % 0 for every base in the locus seq

probe_indcs3 = []; % array where we store index of each probe aligned to target

for i = 1:numel(probe_seqs2_filtered)

    cur_probe = probe_seqs2_filtered{i};
    
    % find where each oligo lands 
    k = strfind(locus_seq, cur_probe);
    
    probe_indcs3 = [probe_indcs3, k];
end

probe_aligned3(probe_indcs3) = 1;

figure
plot(1:length(probe_aligned), probe_aligned);
hold on
plot(1:length(probe_aligned3), probe_aligned3)

%% going to keep only half

% keep_even_idcs = rem(1:numel(probe_seqs2_filtered), 2);
% keep_even_idcs(1:10)= [1]; 
% 
% probe_seqs2_filtered_even = probe_seqs2_filtered(logical(keep_even_idcs), 1);
% 
% numel(probe_seqs2_filtered_even)

%% save
% saveFile = fullfile(saveDir, 'filtered_probes2.csv');

writematrix(cell2mat(probe_seqs2_filtered_even), saveFile);

disp(strcat('Probes were saved to: ', saveFile));

