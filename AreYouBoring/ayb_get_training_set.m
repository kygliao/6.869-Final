% get_training_set 
% Training the data

load VOC2010/person_final.mat;
im = imread('000061.jpg');
disp(size(im));
disp(model);
% TODO: I want to see whether these actually give bounding boxes. but
% imgdetect is not functioning properly at the moment and I'm not sure why.
boxes = imgdetect(im, model, -0.5);
disp(boxes);
keyboard;
[ds, bs] = process(im, model, -0.5);
showboxes(im, [ds,bs]);
% TODO: get the boxes from process
sub_images = ds % this should be the set of detection bounding boxes after clipping
% in ds, each row is a box represented by coordinates [x1 x2 y1 y2]
for i = 1:size(sub_images, 1) %for each box in the set of sub_images
    % TODO: take the bounding box and extract the image from im
    xmin = min(sub_images(i,1), sub_images(i,2));
    ymin = min(sub_images(i,3), sub_images(i,4));
    width = abs(sub_images(i,2) - sub_images(i,1)) + 1;
    height = abs(sub_images(i,4) - sub_images(i,3)) + 1;
    disp(i);
    disp(uint8([xmin ymin width height]));
    I2 = imcrop(im, uint8([xmin ymin width height]));
    
    % TODO: save the im into a training set folder
    % need to test this function
    filename = strcat('ayb_training_set/train', num2str(i));
    filename = strcat(filename, '.jpg');
    imwrite(I2, filename, 'jpg');
end