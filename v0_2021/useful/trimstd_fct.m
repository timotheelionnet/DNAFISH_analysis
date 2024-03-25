function [m,x] = trimstd_fct(data,x)

m = zeros(size(x));

for i=1:numel(x)
    m(i) = trimstd(data,x(i));     
end
end