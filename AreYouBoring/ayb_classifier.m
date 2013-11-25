function [ verdict ] = ayb_classifier(  )
%ayb_classifier takes trained features and an image; outputs "a" or "b"
%   Given training features and an image, the classifier will decide if the
%   image shows a person who is paying attention or bored

addpath(genpath(pwd));

% Initialize variables for calling datasets_feature function
info = load('kai_filelist.mat'); 
datasets = {'demo'}; 
train_lists = {info.train_list};

feature = 'hog2x2'; % let's try gist, hog, and sift; color probably isn't a good idea

c = conf();
c.feature_config.(feature).dictionary_size=20;

% Compute train and test features
% we don't actually need test_labels
test_lists = {info.test_list}; % TODO: replace this with list of images in ayb_tmp
datasets_feature(datasets, train_lists, test_lists, feature, c);

% TODO: fitting it all into an SVM somehow...ugh.
% Trying to interpret fast-additive-svms -- we can likely adapt that code.
% we can temporarily do something funky by scoring each of the images with
% +1 or -1 for attentiveness and then find the sum according to each of the
% features. It's not SVM, but it's a linear combination of stuff. make sure
% we use an odd number of features so that each image has a strictly
% positive or negative score

% Load train and test features
train_features = load_feature(datasets{1}, feature, 'train', c);
test_features = load_feature(datasets{1}, feature, 'test', c);

train_labels = info.train_labels; 
% nearest neighbor computation
[~, nn_idx] = min(sp_dist2(train_features, test_features));
verdict = train_labels(nn_idx);
disp(verdict);

% verdict gives us the matrix of 1xn image categorizations. 
% 1 = attention
% 2 = bored

end

