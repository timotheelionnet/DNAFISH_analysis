function [res,a] = fetch_annotation_from_FlydbID(ID,ID_type)
%Use this function to find the complete annotation from Flydb based on
%either Flybase #, gene symbol, annotation, etc.

%INPUT
%Single request: enter query symbol or ID as a string, then 
%enter the ID_type. ID_type can either be 'Symbol', 'FB#', 'Annotation_ID'
    %examples:
    %res = fetch_annotation_from_FlydbID('FBgn0260799','FB#') ;
    %res = fetch_annotation_from_FlydbID('CG17484','Annotation_ID') ;
    %res = fetch_annotation_from_FlydbID('p120ctn','Symbol') ;
    
%Multiple request: enter query as a cell array of symbols (or IDs, annotations)
%all cell entries must be the same type (e.g. all gene symbols, or all
%FB#,etc). As for single entry, ID_type can either be 'Symbol', 'FB#', 'Annotation_ID'
    %examples:
    %res = res = fetch_annotation_from_FlydbID({'FBgn0023178','FBgn0003068'},'FB#') ;
    %res = fetch_annotation_from_FlydbID({'pdf','per'},'Symbol') ;
    %res = fetch_annotation_from_FlydbID({'CG6496','CG2647'},'Annotation_ID') ;

%OUTPUT
    %res is a cell array with as many rows as queries and 5 columns
    %each row corresponds to the respective query.
    %columns are:
    %1 Gene_Symbol
    %2 Primary_FBgn#
    %3 Secondary_FBgn#
    %4 annotation_ID
    %5 secondary_annotation_ID
    
    %If no match was found, entries are empty.
    
%all matches are case insensitive.    
tic
%% requires the file fbgn_annotation_ID_fb_2016_01.tsv !!!!!!
%downloaded from fly database
%http://flybase.org/static_pages/docs/datafiles.html#accessing_files
%columns as follows:
% ##gene_symbol	primary_FBgn#	secondary_FBgn#(s)	annotation_ID	secondary_annotation_ID(s)

%% Load Flydb annotation file in matlab as a cell array
filename = '/Users/lionnett/Documents/FlyISH/Flydb_annotations/fbgn_annotation_ID_fb_2016_01.tsv';
fid = fopen(filename);
a = textscan(fid,'%s %s %s %s %s','Delimiter','\t','MultipleDelimsAsOne',0);

fclose(fid);
for i=1:numel(a)
    a{i} = a{i}(7:end);
end
%this is the annotation file formatted as a Ngenes x 5 cell array
%Columns are:
%1 Gene_Symbol
%2 Primary_FBgn#
%3 Secondary_FBgn#
%4 annotation_ID
%5 secondary_annotation_ID
a = [a{1},a{2},a{3},a{4},a{5}];


%% check that input has the desired format (string or cell of strings).
if ischar(ID)
    ID = {ID};
end

if ~iscell(ID)
    Symb = [];
    disp('ID needs to be either string or cell array of strings');
    return
end


%% Flydb symbol fetched from URL content
if isempty(ID_type)
    ID_type = 'FB#';
end

switch ID_type 
    case 'FB#'
        refcol = 2;
    case 'Annotation_ID'
        refcol = 4;
    case 'Symbol'
        refcol = 1;
end

%add last coma to fifth column of a() (the column that holds the secondary annotation_ID, looking like "" or
%"CG34226,CG30029,CG18337")
%this ensures that we can loook for the string 'CG30029,'
%if no comma, we would get positive hits for CG30029 and CG300291 for
%instance.

%x: indices of non-zero entries (gene rows that actually have secodnary
%annotation_IDs
%other entries will be ignored to accelerate the search
secondary_idx = cell2mat(cellfun(@size,a(:,5),'UniformOutput',0)); 
secondary_idx = secondary_idx(:,2); 
secondary_idx = find(secondary_idx);
secondary_ID = cellfun(@cat,repmat({2},size(a(secondary_idx,5))),a(secondary_idx,5),repmat({','},size(a(secondary_idx,5))),'UniformOutput',0);

%obsolete
%secondary_ID = regexp(a(:,5),'[\s\,]*','split');

tbIdx = zeros(numel(ID),1);
nohit_Idx = [];
k=0;
for i=1:numel(ID)
    %find the index of the line from table a in which the annotation query was found
    t = find(ismember(a(:,refcol), ID{i})); 

    %if no hit found in primary annotation_ID (column refcol), I add a dummy index (1) as a placeholder (will be deleted
    %later)
    if isempty(t) && strcmp(ID_type,'Annotation_ID')
        
        t = find(~cellfun(@isempty,cellfun(@strfind,secondary_ID,repmat({[ID{i},',']},size(secondary_ID)),'UniformOutput',0)));
        
        
        if isempty(t)
            tbIdx(i) = 1;
            k = k+1;
            nohit_Idx(k) = i;
        else
            %take first result if one secondary annotation points to multiple
            %entries (seems to happen sometimes, not sure this is the best
            %choice to take)
            t = t(1);

            %map index to full table (secondary_ID table has many rows
            %eliminated)
            t = secondary_idx(t);
            tbIdx(i) = t;
        end
    else
        tbIdx(i) = t;
    end
end


res = a(tbIdx,1:5);
res(nohit_Idx,:) = {[]};

t = toc
