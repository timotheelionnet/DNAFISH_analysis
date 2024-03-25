function is_identical = test_whether_2_arrays_are_identical(a1,a2,warning_text)

is_identical = 0;
if isequal(size (a1), size(a2)) || (isvector(a1) && isvector(a2) && numel(a1) == numel(a2))
    if max(abs(diff(a1,a2))) < 1e-8
        is_identical = 1;
    end
end

if ~is_identical && ~isempty(warning_text)
    disp(warning_text);
end
    

end

