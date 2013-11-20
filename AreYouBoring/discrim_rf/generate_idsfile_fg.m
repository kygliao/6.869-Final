function generate_idsfile_fg(trainingpath)

[fid,msg] = fopen(trainingpath);

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


ids_fg = cell(count,1);

[fid,msg] = fopen(trainingpath);

tline = fgetl(fid);
count = 0;
savefile = 'ids_fg.mat';

while ischar(tline)

count = count + 1;


c = tline(1:end-4);

%disp(c);
%disp(c(1:11));
ids_fg{count,1}=c;

tline = fgetl(fid);

end

fprintf('number of training/val fg data is %d\n', count);



fclose(fid);


save(savefile, 'ids_fg');
end
