% create_filelist
% script for creating the proper file list to train and test images

atten_dir = 'set_attention/';
bored_dir = 'set_bored/';

% someone used some interesting .NET extension of MATLAB, but I don't wanna
% go through that. there must be another way to get files...
atten_list = dir(atten_dir);
bored_list = dir(bored_dir);

atten_files = {atten_list.name};
bored_files = {atten_list.name};

[~, a_size] = size(atten_files);
[~, b_size] = size(bored_files);

% This is just here if we want to do some testing. Otherwise, I'm gonna be
% lazy and just collect training data to get training features
%a_x = a_size/2;
%b_x = b_size/2;
%test_labels = ones(1, a_x);
%test_labels(1, a_x + 1 : a_x + b_x) = 2;
%train_labels = test_labels; %if we decide that they should be the same

classes = [1 2];

train_list = [atten_files bored_files];
train_labels = ones(1, a_size);
train_labels(1, a_size + 1 : a_size + b_size) = 2;

%save('traininglist.mat', 'classes', 'test_labels', 'test_list', 'train_labels', 'train_list');
save('traininglist.mat', 'classes', 'train_labels', 'train_list');