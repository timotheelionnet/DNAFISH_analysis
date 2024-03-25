function str2 = parse_and_format_string_list(str)
%input: a list of string separated by comas or smi colons
%output: a cell array of the substrings
%spaces are removed and ignored


%remove spaces
str = strrep(str,' ', '');

%return empty output if empty input
if isempty(str)
    str2 = {[]};
    return
end


%find locations of delimiters
%add artifical delimiters at beginning and end
str = [',',str,','];
x1 = regexp(str,',');
x2 = regexp(str,';');
delimiters_idx = sort(unique([x1,x2]));  

%format strings between delimiters as a cell array of strings
str2 = [];
j=1;
for i=1:numel(delimiters_idx)-1
    tmpstr = str(delimiters_idx(i)+1:delimiters_idx(i+1)-1);
    if ~isempty(tmpstr)
        str2{j} = tmpstr;
        j = j+1;
    end
end
end