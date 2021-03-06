% get_training_set 
% Training the data
LARGE_CLASSROOM_THRESHOLD = 20;
load VOC2010/person_final.mat;
for q = 14:75
    

    IsLargeClassroom=0;

    %load('annotation_using_a_computer.mat');
    file = strcat('ayb_classroom_images/', num2str(q));
    file = strcat(file,'.jpg');
    im = imread(file);

    %disp(size(im));
    %disp(model);
    % TODO: I want to see whether these actually give bounding boxes. but
    % imgdetect is not functioning properly at the moment and I'm not sure why.
    %boxes = imgdetect(im, model, -0.5);
    %disp(boxes);
    [ds, ~] = process(im, model, -0.5);
    %disp(ds);
    showboxes(im, ds);
    % TODO: get the boxes from process
    sub_images = ds; % this should be the set of detection bounding boxes after clipping
    % in ds, each row is a box represented by coordinates [x1 x2 y1 y2]
    numfilters = floor(size(ds, 2)/4);
    for i = numfilters:-1:1 %for each filter
        % TODO: take the bounding box and extract the image from im
        %xmin = min(sub_images(i,1), sub_images(i,2));
        %ymin = min(sub_images(i,3), sub_images(i,4));
        %width = abs(sub_images(i,2) - sub_images(i,1)) + 1;
        %height = abs(sub_images(i,4) - sub_images(i,3)) + 1;
        x1 = ds(:,1+(i-1)*4);
        y1 = ds(:,2+(i-1)*4);
        x2 = ds(:,3+(i-1)*4);
        y2 = ds(:,4+(i-1)*4);
        %disp(i);
        %disp(y2);
        width=x2-x1;
        height= y2-y1;
        %disp(uint8([x1 y1 width height]));
        xmin = uint16(x1);
        ymin = uint16(y1);
        width = uint16(width+1);
        height = uint16(height+1);
        %disp([xmin ymin width height]);
        boxesInfo=[xmin ymin width height];
        %disp(ds-boxesInfo);
        for j = 1:size(ds, 1)
        %I2 = imcrop(im, [x1 y1 width height]);
            %disp(boxesInfo(j,:));
            I2 = imcrop(im,boxesInfo(j,:));
        % TODO: save the im into a training set folder
        % need to test this function
            filename = strcat('ayb_training_set/train',num2str(q));
            filename = strcat(filename,'_');
            filename = strcat(filename, num2str(j));
            filename = strcat(filename, '.jpg');
            imwrite(I2, filename, 'jpg');
        end
        if size(ds, 1)>LARGE_CLASSROOM_THRESHOLD
            IsLargeClassroom=1;

        end
          display('Image:');
          display(q);
          display(IsLargeClassroom);
    end
end