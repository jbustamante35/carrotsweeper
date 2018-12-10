% get the file list of images to trace
FilePath = '/mnt/tetra/nate/IF_gravitropism_examples/';
FileList = {};
FileExt = {'tiff','TIF'};
verbose = 1;
SET = sdig(FilePath,FileList,FileExt,verbose);
%% view set S
s = 1;
close all
M = [];
fM = [];
%% 
para.scales.value = 2;
para.resize.value = .75;
funcIN = @(X)surKur(X,para);
I = double(imread(SET{1}{1}))/255;
R = [2 10];
N = [9 50];
mag = 1;
[BUG,tmpR] = buyEye_NET(I,R,N,mag,false,funcIN,[]);
%%
R = [2 10];
N = [9 50];
mag = 1;
dataXFunc = @(X)buyEye_NET(X,R,N,mag,false,funcIN,[]);
R2 = [10 10];
N2 = [1 50];
mag = 1;
dataYFunc = @(X)buyEye_NET(X,R2,N2,mag,false,[],[]);
[sigX sigY] = generateStemDomain(SET{1},idxL{1},dataXFunc,dataYFunc,R);
%%
inputSize = size(sigX{1},1);
outputSize = 5;
outputMode = 'sequence';
numClasses = 2;
layers = [ ...
    sequenceInputLayer(inputSize)
    lstmLayer(outputSize,'OutputMode',outputMode)
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];

maxEpochs = 5;
options = trainingOptions('sgdm', ...
    'InitialLearnRate',0.001,...
    'MaxEpochs',maxEpochs);


netROW = trainNetwork(sigX(1:10:end),sigY(1:10:end)',layers,options);
%%
tY = netROW.predict(sigX(4000));
%%
I = imread(SET{1}{1});
%[c r V] = impixel(I);
dataXFunc = @(X)buyEye_NET(X,R,N,mag,false,funcIN,netROW);
dataXFunc(I);

%% extract masks
parfor e = 1:numel(SET)
   [idxL{e}] = extractStemPoints(SET{e},false);
end
%% find start points
for e = 1:numel(SET)
    [startPoints{e}] = findStartPoints(idxL{e},SET{e});
end
%% trace the skeletons from each start point to all branch and endpoints
close all
parfor e = 1:numel(SET)
    PATHS{e} = traceSkeletons(SET{e},idxL{e},startPoints{e},false);
end
%% trace curves through time
close all
for e = 1:numel(SET)
    tracePath(SET{e},PATHS{e});
end
%% find start curve
for e = 1:numel(SET)
    initIDX{e} = findInitCurve(SET{e},idxL{e},PATHS{e},startPoints{e});
end
%%
[d] = pathDistance(PATHS{1}{1}{1}{1},PATHS{1}{1}{2}{1});
%%
cnt = 1;
for e = 1:numel(PATHS)
    for sp = 1:numel(PATHS{e})
        for t = 1:numel(PATHS{e}{sp})
            for pth = 1:numel(PATHS{e}{s}{t})
                tmp = PATHS{e}{s}{t}{pth};
                d = diff(tmp,1,1);
                d = sum(d.^2,2).^.5;
                l(cnt) = sum(d);
                cnt =  cnt + 1;
            end
        end
    end
end
%%
for e = 1:numel(SET)
    viewOverlay(SET{e},idxL{e},startPoints{e},initIDX{e},PATHS{e});
end
%%
dM = diff(fM,1,3);
S = mean(abs(dM),3);
close all
%S = std(udM,1,3);
imshow(S,[])
S1 = mean(S,1);
figure;plot(S1);

