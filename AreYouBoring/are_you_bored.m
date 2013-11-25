function [ output_args ] = AreYouBored( image )
%AreYouBored provides the user with information about the audience
%according to an image taken during the presentation
%   input arguments:
%       image file name
addpath(genpath(pwd));

load VOC2010/person_final.mat;

% initialize the count of number people paying attention and number bored
attention = 0;
bored = 0; 

im = imread(image);
[ds, ~] = process(im, model, -0.5);

% code taken from "ayb_get_training_set.m"
LARGE_CLASSROOM_THRESHOLD = 20;
IsLargeClassroom=0;
sub_images = ds; % this should be the set of detection bounding boxes after clipping
% in ds, each row is a box represented by coordinates [x1 x2 y1 y2]
numfilters = floor(size(ds, 2)/4);
for i = numfilters:-1:1 
    x1 = ds(:,1+(i-1)*4);
    y1 = ds(:,2+(i-1)*4);
    x2 = ds(:,3+(i-1)*4);
    y2 = ds(:,4+(i-1)*4);
    width=x2-x1;
    height= y2-y1;
    xmin = uint16(x1);
    ymin = uint16(y1);
    width = uint16(width+1);
    height = uint16(height+1);
    boxesInfo=[xmin ymin width height];
    for j = 1:size(ds, 1)
        I2 = imcrop(im,boxesInfo(j,:));
        filename = strcat('ayb_tmp/sub_',num2str(j));
        filename = strcat(filename, '.jpg');
        imwrite(I2, filename, 'jpg');
    end
    if size(ds, 1)>LARGE_CLASSROOM_THRESHOLD
        IsLargeClassroom=1;

    end
      
end

display(IsLargeClassroom);

verdict = ayb_classifier();
for i=1:size(verdict, 2)
    if verdict(1,i) == 1
        attention = attention + 1
    else
        if verdict(1,i) == 2
            bored = bored + 1
        end
    end

end

% not letting me commit again. =(

%TODO: return the actions and interpretation for the user
end

