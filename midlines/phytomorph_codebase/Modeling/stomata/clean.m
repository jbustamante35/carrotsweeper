D = readtext('/home/nate/Downloads/SD_2016_genotypeNames.csv');
%D = readtext('/home/nate/Downloads/sd_2017_names.csv');
%%
FilePath = '/mnt/tetra/nate/stomataTopoData/Accessions_2016/';
FileList = {};
FileExt = {'nms'};
tic
FileList = gdig(FilePath,FileList,FileExt,1);
toc
%%
X = table;
cnt = 1;
for e = 1:1000%numel(FileList)
    [pth,nm,ext] = fileparts(FileList{e});
    ridx = find(strcmp(D(:,2),strrep(['16' nm],'','')))
    if ~isempty(ridx)
        X(cnt,'filename') = {FileList{e}};
        X(cnt,'Count') = {D{ridx,3}};
        Y(cnt) = D{ridx,3};
        cnt = cnt + 1;
    end
end
%%
I = imread(FileList{1});
layers = [imageInputLayer(size(I),'Normalization','none');
          convolution2dLayer(51,3);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          convolution2dLayer(21,3);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          convolution2dLayer(3,3);
          reluLayer();
          maxPooling2dLayer(2,'Stride',2);
          fullyConnectedLayer(1)
          regressionLayer];
options = trainingOptions('sgdm','InitialLearnRate',0.0001,'Plots','training-progress','MaxEpochs',45,'ExecutionEnvironment','parallel');
trainedNet = trainNetwork(X,layers,options);



%%  remove Q1
for e = 1:size(D,1)
    if contains(D{e,11},'Q1')
        rm(e) = true;
    else
        rm(e) = false;
    end
end
D(rm,:) = [];
%%
%genoData = D(2:end,end);
genoData = D(2:end,10);
genoNames = unique(genoData);
%phenoData = cell2mat(D(2:end,4));
phenoData = cell2mat(D(2:end,7));
sampleData = D(2:end,4);
%%
ksdensity(phenoData)
%%
close all
gm = {};
AIC = [];
OUTPUT = {};
for e = 1:numel(genoNames)
    fidx = find(strcmp(genoData,genoNames{e}));
    subP = phenoData(fidx,:);
    options = statset('Display','off');
    samp = sampleData(fidx);
    
    for g = 1:2
        flag = 0;
        cnt = 1;
        while flag == 0
            try
                %isoutlier(subP,'method','gesd')
                
                gm{e,g} = fitgmdist(subP,g,'Options',options,'Replicates',10,'RegularizationValue',.01);
                AIC(e,g) = gm{e,g}.AIC;
                %AIC(e,g) = gm{e,g}.BIC;
                %AIC(e,g) = -gm{e,g}.NegativeLogLikelihood;
                flag = 1;
            catch
                flag = 0;
                cnt = cnt + 1;
                if cnt > 20
                    flag = 1;
                    AIC(e,g) = inf;
                end
                
            end
        end
        
    end
    
    
    [Y,X] = ksdensity(subP);
    
    %plot(X,Y,'b');
    %hold on
    [~,midx] = min(AIC(e,:));
    A = gm{e,midx};
    kidx = A.cluster(subP);
   
    [~,MIDX] = max(A.ComponentProportion);
    k = kidx == MIDX;
    
    CL = {'k'};
    for s = 1:numel(A.ComponentProportion)
        Yi = normpdf(X,A.mu(s),A.Sigma(s).^.5);
     %   plot(X,Yi,CL{1});
        if s == MIDX
       %    plot(X,Yi,'r'); 
            valueP = A.mu(s);
        end
    end
    
    %title([genoNames{e} '-->' num2str(e)])
    %pause(.01);
    %waitforbuttonpress
    
    %close all
    %{
    ksdensity()
    drawnow
    hold on
    %}
    
    OUTPUT{e,1} = genoNames{e};
    OUTPUT{e,2} = valueP;
    OUTPUT{e,3} = mean(subP);
    e
    numel(genoNames)
end
%%
close all
good = find(AIC(:,1) < AIC(:,2));
good = min(AIC,[],2);
for e = 1:numel(good)
    fidx = find(strcmp(genoData,genoNames{good(e)}));
    subP = phenoData(fidx,:);
   
    
    ksdensity(subP)
    drawnow
    hold on
    waitforbuttonpress
end
