function [ output_args ] = AreYouBored( image )
%AreYouBored provides the user with information about the audience
%according to an image taken during the presentation
%   input arguments:
%       image file name

load VOC2010/person_final.mat;
im = imread(image);
bbox = process(im, model, -0.5);
pd = bbox(1); %bounding boxes for people detected

% TODO: create a hashmap of the actions we care about to # people
% look at "containers.Map"
for i = 1:size(pd, 1)
    % TODO: Get image from boundaries
    xmin = pd(i, 1);
    ymin = pd(i, 3);
    width = pd(i, 2) - pd(i, 1) + 1;
    height = pd(i, 4) - pd(i,3) + 1;
    sub_image = imcrop(im, [xmin, ymin, width, height]);
    % TODO: run image categorization according to the sub_image
    
    % once categorization is complete, update the hashmap
end

%TODO: return the actions and interpretation for the user
end

