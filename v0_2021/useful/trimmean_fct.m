function [m,x] = trimmean_fct(data,x)

m = zeros(size(x));

for i=1:numel(x)
    m(i) = trimmean(data,x(i));     
end
end







