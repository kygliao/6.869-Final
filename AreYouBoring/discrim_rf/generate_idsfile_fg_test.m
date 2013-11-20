function generate_idsfile_fg_test(testpath)

[fid,msg] = fopen(testpath);

tline = fgetl(fid);
str = 'randomstring';
count = 0;
result = [];
hi ='randomforest/llc_extraction/dataset_fg';
name = 'VOC Action Classification fg';

while ischar(tline)

count = count + 1;

tline = fgetl(fid);

end
fclose(fid);


ids_fg_test = cell(count,1);
savefile = 'ids_fg_test.mat';

[fid,msg] = fopen(testpath);

tline = fgetl(fid);
count = 0;

while ischar(tline)

count = count + 1;

c = tline(1:end-4);

%disp(c);
ids_fg_test{count,1}=c;

tline = fgetl(fid);

end

fprintf('number of test fg data is %d\n', count);


fclose(fid);


save(savefile, 'ids_fg_test');
end
