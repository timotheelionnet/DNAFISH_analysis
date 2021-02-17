function m = trimstd(data,x)

data = sort(data,'ascend');

i1 = round((100 - x)/200*numel(data));
i2 = round((100 + x)/200*numel(data));

data = data(i1:i2);

m = std(data(:));

end

