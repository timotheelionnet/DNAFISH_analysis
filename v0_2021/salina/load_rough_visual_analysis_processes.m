%load_rough_visual_analysis_processes

data_xls = '/Users/lionnett/Documents/Salina Long/PDF TS vs ZT/initial_conversion_log_parsed.xlsx';
[num,txt,raw] = xlsread(data_xls);
%data format: excel spreadsheet with in sheet 1; 
%one header line.
%col 1: (not used) minimum intensity across entire 32bit stack 
%col 2: (not used) maximum intensity across entire 32bit stack
%col 3: file name
%col 4: number of channels
%col 5: number of neurons visible in 2D projection
%col 6: comments on processes from 2D projection - contain parsable info!
   %example format is thick-47;thin-46;blob;
   %this means one 47pix long thick process, one 46pix long thin process,
   %one blob seen on 2D projection.
   %see process_terms and other_terms for the vocabulary used.
   
%col 7: (not used in this early version) number of neurons in stack (from 3D dataset inspection)
%col 8: (not used in this early version) number of processes traced in 3D
%col 9: comments from 3D tracing
%col 10: some stats on each condition


%% vocabulary used in column 6 to generate parsable comments on neuron processes.
process_terms = {'thin','thick','faint','thinbranched','sparse'};
other_terms = {'diffuse_spots','blob','no','isolated_spots'};

%% extract data in usable format
%ZT time of ith file
zt = zeros(size(num,1),1); 

%1 if ith file contains sLNvs only
sLNv = false(size(num,1),1);

%1 if ith file contains lLNvs only
lLNv = false(size(num,1),1);

%number of (2D) visible neurons in file i
n_visible_neurons = zeros(size(num,1),1);

%number of (2D) visible processes in file i
n_visible_processes = zeros(size(num,1),1);

%array of processes lengths in file i
process_length = [];

%number of "other" comments in file ith (mainly for checking purposes).
n_others = zeros(size(num,1),1);

for i=1:numel(zt)
   %extract ZT from file name (col 3)
   %line index takes in account header line
   tmpstr = raw{i+1,3};
   ks = strfind(tmpstr,'/zt');
   tmpstr2 = tmpstr(ks+3:end);
   ks = strfind(tmpstr2,'_');
   tmpstr3 = tmpstr2(1:ks(1)-1);
   zt(i) = str2double( tmpstr3 ); 
   
   %extract neuron type from file name (col 3) 
   %line index takes in account header line
   sLNv(i) = ~isempty(strfind(tmpstr,'sLNv')) || ~isempty(strfind(tmpstr,'sLNV'));
   lLNv(i) = ~isempty(strfind(tmpstr,'lLNv')) || ~isempty(strfind(tmpstr,'lLNV'));
   tmpstr = raw{i+1,6};
   tmpstr = strsplit(tmpstr,';');
   
   %collect number of (2D) visible neurons
   n_visible_neurons(i) = raw{i+1,5};
   
   %compile approximate process lengths from comments (col 6)
   %example format is thick-47;thin-46;blob;
   %this means one 47pix long thick process, one 46pix long thin process,
   %one blob seen on 2D projection.
   process_length{i} = [];
   for j=1:numel(tmpstr)
       for k=1:numel(process_terms)
        if ~isempty(strfind(tmpstr{j},process_terms{k}))
            n_visible_processes(i) = n_visible_processes(i)+1;
            l = strfind(tmpstr{j},'-');
            
            process_length{i} = [process_length{i},str2double(tmpstr{j}(l+1:end))];
        end
       end
       for k=1:numel(other_terms)
        if ~isempty(strfind(tmpstr{j},other_terms{k}))
            n_others(i) = n_others(i)+1;
        end
       end
   end
   disp(['i = ',num2str(i),'; found ',num2str(n_visible_processes(i) + n_others(i)),' out of ',num2str(numel(tmpstr)-1),' annotations']);
end

%% compile numbers
%sorted array of ZT conditions
zt_c = sort(unique(zt));
%self explanatory metrics from 2D observations ranked according to zt_c
%order
%avg_nproc_per_neuron
%avg_nproc_per_sLNv
%avg_nproc_per_lLNv
%avg_proc_length

for i=1:numel(zt_c)
    idx = zt == zt_c(i);
    avg_nproc_per_neuron(i) = sum(n_visible_processes(idx))/sum(n_visible_neurons(idx));
    avg_nproc_per_sLNv(i) = sum(n_visible_processes(idx & sLNv))/sum(n_visible_neurons(idx & sLNv));
    avg_nproc_per_lLNv(i) = sum(n_visible_processes(idx & lLNv))/sum(n_visible_neurons(idx & lLNv));
    
    x = process_length(idx);
    p = [];
    for j=1:numel(x)
        p = [p,process_length{j}];
    end
    avg_proc_length(i) = mean(p);
   
end

%plot data around the clock
figure('Name','Average Number of Processes per Neuron'); 
hold;
plot(zt_c,avg_nproc_per_neuron,'DisplayName','all neurons');
plot(zt_c,avg_nproc_per_sLNv,'DisplayName','sLNvs');
plot(zt_c,avg_nproc_per_lLNv,'DisplayName','lLNv');

figure('Name','Average Process Length');
hold;
plot(zt_c,avg_proc_length);


