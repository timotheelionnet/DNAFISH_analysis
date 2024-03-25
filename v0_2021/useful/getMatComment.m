function txt=getMatComment(x)
fid=fopen(x);
txt=char(fread(fid,[1,140],'*char'));
txt=[txt,0];
txt=txt(1:find(txt==0,1,'first')-1);
end