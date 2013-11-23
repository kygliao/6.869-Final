% generates data on two annular rings 
% linear classifier will fail to train well
if(~exist('pwl_sgd.mexa64','file') || ~exist('pwl_sgd.mex','file')),
    fprintf('compiling PWLSGD...\n');
    mex pwl_sgd.cc;
end
mu1    = [0  0];
sigma1 = [1 0; 0 1];
sigma2 = [2 0; 0 2];
npts   = 5000;

%% Training images
%% @TODO figure out how to load images
%% Strategy: divide images into two folders.
%% If the person is paying attention, put the image in one folder
%% else, put it in the other.
%% This will determine if the data point is a 1 or -1.

% How to train TODO
% 1) take in the image, and extract the features
% 2) pass the features into the algorithm to be trained
% X = features at nth point (some matrix)
% We probably don't need p1 and p2 since we don't need those for T
% TODO: call code that extracts features
% 
X      = mvnrnd(mu1,sigma1,npts)';
p1     = mvnpdf(X',mu1,sigma1);
p2     = mvnpdf(X',mu1,sigma2);
T      = zeros(1,npts);
keyboard;
% for loop to iterate through files
% if file is from a specific folder, assign it the corresponding number
% TODO: syntax for using T? set position n to equal -1 or 1 T(n) = 1 or -1;
% if check filename == pos
%   T(n) = 1;
% else
%   T(n) = -1;
% end
T(p1 >= p2) = 1;
T(p1 < p2) = -1;
%% testing images
%% @TODO figure out how to load images
%% Same as above.
X_      = mvnrnd(mu1,sigma1,npts)';
p1      = mvnpdf(X_',mu1,sigma1);
p2      = mvnpdf(X_',mu1,sigma2);
T_      = zeros(1,npts);
T_(p1 >= p2) = 1;
T_(p1 < p2) = -1;


%% plot the data

X  = single(X) ; T  = single(T);
X_ = single(X_); T_ = single(T_);

%learn classifier using stoc gradient descent
lambda         = 0.0001;
k              = 10;
niters         = 10000;
class_1_weight = 1;
param.ndivs    = 20;
param.use_pwl  = 1;
init_iter      = 0;
use_sqrt       = 1;


[Xe,seg_data] = encode_data(X,param);
Xe_           = encode_data(X_,param,seg_data);
tic
w  = pwl_sgd(X,seg_data,T,k,...
             niters,lambda,param.use_pwl,init_iter,...
             use_sqrt,class_1_weight);
toc
pT =  w'*Xe_;  % predicted scores on test data

figure; 
hold on; axis equal;box on;
title('Training Data')
plot(X(1,T > 0), X(2,T > 0),'+g');
plot(X(1,T < 0), X(2,T < 0),'+r');

figure;
hold on; axis equal;box on;
plot(X_(1,pT > 0), X_(2,pT > 0),'+g');
plot(X_(1,pT < 0), X_(2,pT < 0),'+r');
%find misclassified ones
ii = find(sign(pT).*sign(T_) < 0);
plot(X_(1,ii), X_(2,ii),'ok');
error_rate = length(ii)*100/length(pT);
title(sprintf('Classification on Test Data using PWLSGD (%.2f%% error)',error_rate));

%% optionally train a liblinear model

% add path to liblinear
doLIBLINEAR = 1;
if(doLIBLINEAR)
    fprintf('\nTraining LIBLINEAR models on encoded data\n');
    tic
    model = train(double(T'),Xe,'-s 3 -c 10','col');
    toc;
    [label,acc,dec_values] = predict(double(T_'),Xe_',model);

    fprintf('\n\nComparing Test Accuracy:\n\t          PWLSGD : %.2f %% \n\tphi2 + LIBLINEAR : %.2f %%\n',100-error_rate,acc);
end







