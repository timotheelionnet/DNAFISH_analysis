function ts_corr = correct_for_bleachdecay(ts,decay,time)
npts = size(ts,1);
ts_corr = ts;
for i=1:npts
    ti = find(time == ts(i,4));
    ts_corr(i,3) = ts(i,3)/decay(ti);
end

end
