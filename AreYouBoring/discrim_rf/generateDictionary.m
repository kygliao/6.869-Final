function [dictionary] = generateDictionary(data, config)

numImages = config.dictionary.numImages;
numDescriptors = config.dictionary.numDescriptors;
dictionarySize = config.dictionary.size;

imageIdx = randi(length(data.tr), numImages, 1);
patchIdx = zeros(numImages, 1);
siftCount = zeros(numImages, 1);
flipIdx = ones(numImages, 1);
if(config.includeFlippedImages)
    flipIdx = randi(2, numImages, 1);
end

for i=1:numImages
  disp(['generateDictionary: ' num2str(i) ' of ' num2str(numImages)]);
  currImage = data.tr(imageIdx(i)).annotation;
  inputFile = fullfile(config.outputFolder, data.name, 'sift', currImage.folder, [currImage.filename(1:end-4) config.flipSuffix{flipIdx(i)} '.mat']);
  tempData = load(inputFile);

  patchIdx(i) = randi(length(tempData.patchSizes), 1);
  siftCount(i) = length(tempData.data{patchIdx(i)});
end

siftDescriptors = zeros(sum(siftCount), size(tempData.data{1}, 2));
currentIdx = 1;

for i=1:numImages  
  currImage = data.tr(imageIdx(i)).annotation;
  inputFile = fullfile(config.outputFolder, data.name, 'sift', currImage.folder, [currImage.filename(1:end-4) config.flipSuffix{flipIdx(i)} '.mat']);
  tempData = load(inputFile);
  siftDescriptors(currentIdx:currentIdx+siftCount(i)-1, :) = tempData.data{patchIdx(i)};
  currentIdx = currentIdx + siftCount(i);
end

ndata = size(siftDescriptors, 1);
if(ndata > numDescriptors)
  siftDescriptors = siftDescriptors(randsample(ndata, numDescriptors), :);
end

[dictionary, ~] = litekmeans(siftDescriptors', dictionarySize);

dictionary = dictionary';
save([config.outputFolder '/' data.name '/dictionary_' num2str(dictionarySize) '.mat'], 'dictionary');