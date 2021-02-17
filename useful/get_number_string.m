function istr = get_number_string(i,ndigits)
%builds a string from a number i adding zeros before the number to reach a
%total size of ndigits.
%example:
%get_number_string(3,4) yields the string '0003'

if i ==0
    istr = repmat('0',1,ndigits);
    return;
end


nzeros = ndigits - floor(log(i)/log(10)) - 1;
%to correct a rounding error from matlab that sometimes results in:
%floor(log(1000)/log(10)) = 2

if i == 1000
    nzeros = ndigits -4;
end

if nzeros > ndigits
    disp('number has more digits than expected!');
    istr = num2str(i);
    return;
end

if nzeros == 0
    istr = num2str(i);
else
    istr = [repmat('0',1,nzeros),num2str(i)];
end

clear('nzeros');
end