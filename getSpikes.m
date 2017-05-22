function [kullbackError, svmError] = getSpikes()

warning('off','all');

v1List = cell(8,1);
mtList = cell(8,1);
ori = 0;
CreateMethod = 'symmetricAnisotropic';
for i = 1:length(v1List) 
    v1List{i}.center = createV1Cells('center', ori);
    v1List{i}.surround = createV1Cells('surround', ori);
    
    mtList{i}.center = createMTCells('center', ori, CreateMethod);
    mtList{i}.surround = createMTCells('surround', ori, CreateMethod);
    ori = ori + pi/4;
end 

[training, test] = createSamplesForTraining('C:\Users\Sevim Caliskan\Desktop\vision information processing\DVS128 Records\out\videos\', v1List, mtList, 3, [80 80]);
labels = cell(length(test),1);
for i = 1: length(labels)
    labels{i} = test(i).label;
end

kullback = applyKullbackLeiblerDivergence(test, training);
kullbackError = zeros(length(test),1);
for i = 1:length(kullbackError)
    if strcmp(kullback{i}, test(i).label)==1
        kullbackError(i) = 1;
    end
end
kullbackError = 1 - sum(kullbackError)/length(kullbackError);

svmResults = applySVM(test, training);
svmR = cell(length(svmResults)/210,1);
for i = 1:length(svmR)
    interval = ((i-1)*210+1):(i*210);
    temp = svmResults(interval);
    numberof0s = sum(temp==0);
    numberof1s = sum(temp==1);
    numberof2s = sum(temp==2);
    if(numberof0s>numberof1s && numberof0s>numberof2s)
        svmR{i} = 'paper';
    elseif(numberof1s>numberof0s && numberof1s>numberof2s)
        svmR{i} = 'scissor';
    elseif(numberof2s>numberof0s && numberof2s>numberof1s)
        svmR{i} = 'scissor';
    end
end

svmError = zeros(length(test),1);
for i = 1:length(svmError)
    if strcmp(svmR{i}, test(i).label)==1
        svmError(i) = 1;
    end
end
svmError = 1 - sum(svmError)/length(svmError);

end

function spikeList = getMTResponsesMean(frames, mtList, duration, centerSize)
spikeList = [];
for i = 1:length(mtList)
    spikes = get1MTResponseWithV1s(frames, mtList{i}, duration, centerSize);
    spikeList = cat(3, spikeList, spikes);
end
spikeList = sum(spikeList,3);
spikeList = (spikeList - min(min(spikeList)))./(max(max(spikeList)) - min(min(spikeList)));
end

function responses = getV1ResponsesCollection(frames, v1List, duration, centerSize)
    responses = cell(length(v1List), 1);
    parfor i = 1:length(v1List)
        v1Center = v1List{i}.center;
        v1Surround = v1List{i}.surround;
        response = combineV1CenterSurround(frames, v1Center, v1Surround, duration, centerSize);
        responses{i} = response;
    end
end

function spikes = get1MTResponseWithV1s(responses, mt, duration, centerSize)
    spikes = [];
    ori = 0;
    for i = 1:length(responses)
        frames = responses{i};
        response = combineMTCenterSurround(frames, mt.center, mt.surround, duration, ori, centerSize);
        ori = ori + pi/4;
        spikes = cat(3, spikes, response);
    end
    spikes = sum(spikes,3);
    spikes = (spikes - min(min(spikes)))./(max(max(spikes)) - min(min(spikes)));
end

function spikes = combineV1CenterSurround(frames, v1Center, v1Surround, duration, centerSize)
    s = size(frames);
    x_length = s(1);
    y_length = s(2);
    center_x = ceil(x_length/2);
    center_y = ceil(y_length/2);
    
    xCoords = center_x - floor(centerSize(1)/2):1: center_x + floor(centerSize(1)/2);
    yCoords = center_y - floor(centerSize(2)/2):1: center_y + floor(centerSize(2)/2);
    centerResponse = getV1SpikesForOneSample(frames, v1Center, duration, xCoords, yCoords);
    surroundResponse = getV1SpikesForOneSample(frames, v1Surround, duration, 1:size(frames,1), 1:size(frames,2));
    surroundResponse(xCoords, yCoords) = centerResponse;
    spikes = surroundResponse;
    
end


function spikes = combineMTCenterSurround(frames, mtCenter, mtSurround, duration, v1Selectivity, centerSize)
    s = size(frames);
    x_length = s(1);
    y_length = s(2);
    center_x = ceil(x_length/2);
    center_y = ceil(y_length/2);
    
    xCoords = center_x - floor(centerSize(1)/2):1: center_x + floor(centerSize(1)/2);
    yCoords = center_y - floor(centerSize(2)/2):1: center_y + floor(centerSize(2)/2);
    centerResponse = getMTSpikesForOneSample(frames, mtCenter, v1Selectivity, duration, xCoords, yCoords);
    surroundResponse = getMTSpikesForOneSample(frames, mtSurround, v1Selectivity, duration, 1:size(frames,1), 1:size(frames,2));
    surroundResponse(xCoords, yCoords) = centerResponse;
    spikes = surroundResponse;
    
end

function v1 = createV1Cells(place, orientation)
   if strcmp('center', place)==1
      v1 = V1ComplexCell([3 3], 1, orientation, 0.3);
   else
      v1 = V1ComplexCell([7 7], 1, orientation, 0.3);
   end
end

function mt = createMTCells(place, orientation, CreateMethod)
   if strcmp('center', place)==1
      mt = MTCell(orientation, 0.8, CreateMethod, 7, 0);
   else
      mt = MTCell(orientation, 0.4, CreateMethod, 15, 0);
   end
end

function spikes = getV1SpikesForOneSample(frames, v1, duration, xCoords, yCoords)
    spikes = v1.Integrate(frames, duration);
    spikes = spikes(xCoords, yCoords);
end


function spikes = getMTSpikesForOneSample(frames, mt, v1Selectivity, duration, xCoords, yCoords)
    mt.setKernel(v1Selectivity);
    spikes = mt.Integrate(frames, duration);
    spikes = spikes(xCoords, yCoords);
end


function [training, test] = createSamplesForTraining(folder, v1List, mtList, duration, centerSize)
files = dir(strcat(folder, '*.avi'));
cd(folder)
for i = 1:length(files)
    if ~isempty(strfind(files(i).name, 'paper'))
        files(i).label = 'paper';
    elseif ~isempty(strfind(files(i).name, 'scissor'))
        files(i).label = 'scissor';
    elseif ~isempty(strfind(files(i).name, 'rock'))
        files(i).label = 'rock';
    end
end

training = struct('name', {}, 'data', {}, 'label', {});
test = struct('name', {}, 'data', {}, 'label', {});
trainingLength = round(length(files)/4);
trainingIndices = randperm(length(files));
trainingIndices = trainingIndices(1:trainingLength);
for i = 1:length(files)
    frames = videoToFrames(files(i).name);
    cd('C:\Users\Sevim Caliskan\Desktop\vision information processing')
    responses = getV1ResponsesCollection(frames, v1List, duration, centerSize);
    data = getMTResponsesMean(responses, mtList,duration, centerSize);
    if(isempty(find(trainingIndices==i)))
        trainingElement = struct('name', {files(i).name}, 'data', {data}, 'label', {files(i).label});
%         trainingElement.name = files(i).name;
%         trainingElement.data = data;
%         trainingElement.label = files(i).label;
        training = cat(1, training, trainingElement);
    else
        testElement = struct('name', {files(i).name}, 'data', {data}, 'label', {files(i).label});
%         testElement.name = files(i).name;
%         testElement.data = data;
%         testElement.label = files(i).label;
        test = cat(1, test, testElement);
    end
    cd(folder)
end

cd('C:\Users\Sevim Caliskan\Desktop\vision information processing')

end


function triangleResults = applyKullbackLeiblerDivergence(test, training)
   triangleResults = cell(length(test),1);
   for i = 1:length(test)
       minDistance = Inf;
       for j = 1:length(training)
           distance = getTotalEntropyDistance(test(i).data, training(j).data);
           if minDistance>distance
               minDistance = distance;
               triangleResults{i} = training(j).label;
           end
       end
   end
   
end

function svmR = applySVM(test, training)
   x = [];
   t = [];
   for i = 1:length(training)
       features = training(i).data;
       x = [x; features];
       if(strcmp(training(i).label, 'paper')==1)
           temp = zeros(210, 1); temp(:) = 0;
           t = [t;temp];
       elseif(strcmp(training(i).label, 'scissor')==1)
           temp = zeros(210, 1); temp(:) = 1;
           t = [t;temp];
       elseif(strcmp(training(i).label, 'rock')==1)
           temp = zeros(210, 1); temp(:) = 2;
           t = [t;temp];
       end
   end
   Mdl = fitcecoc(x,t);
   
   
   
   x = [];
   t = [];
   for i = 1:length(test)
       features = test(i).data;
       x = [x; features];
       if(strcmp(test(i).label, 'paper')==1)
           temp = zeros(210, 1); temp(:) = 0;
           t = [t;temp];
       elseif(strcmp(test(i).label, 'scissor')==1)
           temp = zeros(210, 1); temp(:) = 1;
           t = [t;temp];
       elseif(strcmp(test(i).label, 'rock')==1)
           temp = zeros(210, 1); temp(:) = 2;
           t = [t;temp];
       end
   end
   
   svmR = predict(Mdl, x);
end

function diff = getTotalEntropyDistance(img, ref)
    imgHist = calcEntropyFromHist(img);
    refHist = calcEntropyFromHist(ref);
    diff = abs(imgHist-refHist);
   
end

function entropy = calcEntropyFromHist(input)
    [h, x] = hist(input,min([max([0.1*size(input,1)*size(input,2),10]),50]));

    zero = find(h == 0);
    for k = 1:length(zero)
        h(zero(k)) = 1;
    end

    h = h/sum(h);
    step = x(2)-x(1);
    entropy = log(step)-sum(h.*log(h));
end

function frames = videoToFrames(pathToFile)
v = VideoReader(pathToFile);
n = v.NumberOfFrames;
frames = [];
for i = 1:n
    frame = read(v, i);
    frame = rgb2gray(frame);
    frame = imresize(frame, [210 210]);
    frames = cat(3, frames, frame);
end
end