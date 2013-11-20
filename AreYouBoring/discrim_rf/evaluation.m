function evaluation
addpath('VOCcode');

[result, classVec1, classVec2] = dttest();
numImages = size(result, 1);
numClasses = size(result, 2);
numTrees = size(result, 3);
numTrees
numClasses
for i = 1 : numTrees
    for j = 1 : numClasses
        [rec, prec, ap(j)] = ProcessClass(result(:, :, i), classVec1, classVec2, j);
    end
    %disp(['The overall mAP until the ' num2str(i) '-th tree: ' num2str(mean(ap(1:10)))]);
    disp(['The overall mAP until the ' num2str(i) '-th tree: ' num2str(mean(ap(1:numClasses)))]);
end

save('../CCodeV16/train_confidence_rf.mat', 'result');
save('../CCodeV16/train_classVec1.mat', 'classVec1');
save('../CCodeV16/train_classVec2.mat', 'classVec2');
for i = 1 : numClasses
    disp(['  The AP for the ' num2str(i) '-th class is: ' num2str(ap(i))]);
end



function [rec, prec, ap] = ProcessClass(result, classVec1, classVec2, currClass)

numImages = size(result, 1);
gt = ones(numImages, 1) .* -1;
gt(find(classVec1 == currClass)) = 1;
gt(find(classVec2 == currClass)) = 1;
out = result(:, currClass);

[so, si] = sort(-out);
tp = gt(si) > 0;
fp = gt(si) < 0;

fp = cumsum(fp);
tp = cumsum(tp);
rec = tp / sum(gt > 0);
prec = tp ./ (fp + tp);
ap = VOCap(rec, prec);
