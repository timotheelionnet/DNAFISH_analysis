function register_mp_and_stack(dir)


i1 = 'mp3.tif';
i2 = 'mp35.tif';
i3 = 'mp5.tif';

f1 = 'mp3reg.tif';
f2 = 'mp35reg.tif';
f3 = 'mp5reg.tif';

nst1 = 'st3.tif';
nst2 = 'st35.tif';
nst3 = 'st5.tif';

nst1f = 'st3reg.tif';
nst2f = 'st35reg.tif';
nst3f = 'st5reg.tif';


mp1 = timtiffread([dir,i1]);
mp2 = timtiffread([dir,i2]);
mp3 = timtiffread([dir,i3]);

[mp1f,mp2f,mp3f,d12,d23,d13,indices] = register_3images(mp1,mp2,mp3,10);

save_as_tiff(mp1f,[dir,f1]);
save_as_tiff(mp2f,[dir,f2]);
save_as_tiff(mp3f,[dir,f3]);

st1 = timtiffread([dir,nst1]);
st1f = st1(indices(1,1,1):indices(2,1,1),indices(1,2,1):indices(2,2,1),:);
save_as_tiff(st1f,[dir,nst1f]);
clear('st1','st1f');

st2 = timtiffread([dir,nst2]);
st2f = st2(indices(1,1,2):indices(2,1,2),indices(1,2,2):indices(2,2,2),:);
save_as_tiff(st2f,[dir,nst2f]);
clear('st2','st2f');

st3 = timtiffread([dir,nst3]);
st3f = st3(indices(1,1,3):indices(2,1,3),indices(1,2,3):indices(2,2,3),:);
save_as_tiff(st3f,[dir,nst3f]);
clear('st3','st3f');

clear('mp1','mp2','mp3','mp1f','mp2f','mp3f','d12','d23','d13','indices');


end