function [P] = getImagePatch(I,oPath,disp)
    % example of square patch
    [x1 x2] = ndgrid(linspace(-30,30,100),linspace(-30,30,100));
    % disk
    rad = 30;
    numrad = 50;
    [R T] = ndgrid(linspace(0,rad,numrad),linspace(-pi,pi,100));
    X = R.*cos(T);
    Y = R.*sin(T);
    %
    skip = 1;
    toKeep = 15;
    close all
    
    %
    ri = (rad+1):skip:(size(I,1)-rad);
    ci = (rad+1):skip:(size(I,2)-rad);
    P = zeros(numel(ri),numel(ci),toKeep*numrad);
    indexR = 1;
    indexC = 1;
    
    for r = ri
        indexC = 1;
        for c = ci
            %{
            tic
            v = interp2(I,X+r,Y+c);
            toc*numel(I)/60
            %}
            
            tic
            v = ba_interp2(I,X+r,Y+c);
            v = bsxfun(@minus,v,mean(v,2));
            v = fft(v,[],2);
            v = abs(v);
            v = v(:,1:toKeep);
            P(indexR,indexC,:) = v(:);
            toc*numel(I)/60
            
            
            if disp
                imshow(I,[]);
                hold on
                plot(Y(:)+c,X(:)+r,'.');
                imshow(bindVec(v),[])
                hold off
                pause(.1)
                drawnow
                %mesh(v)
            end
            indexC = indexC + 1;
        end
        indexR = indexR + 1;
    end
end


%{

    dataPath = ['/iplant/home/leakey_cyverse/maizeData/stomataTopoData/RILs%'];
    CMD = ['iquest --no-page "select COLL_NAME,DATA_NAME where COLL_NAME like ''' dataPath '''"']; 
    [o,r] = system(CMD);
    [r] = parseRecords(r);
    FileExt = {'nms'};
    nmsFileList = {};
    for e = 1:numel(r)
        [p,nm,ext] = fileparts(r(e).DATA_NAME);
        if any(strcmp(ext(2:end),FileExt))
            nmsFileList{end+1} = [r(e).COLL_NAME filesep r(e).DATA_NAME];
        end
    end


    fileName = ['/iplant/home/leakey_cyverse/maizeData/stomataTopoData/RILs/624 leaf4-4.nms'];
    P = zeros([452 452 750 30]);
    parfor e = 1:30
        I = imread(nmsFileList{e});
        P(:,:,:,e) = getImagePatch(I,'',false);
    end

    osz = size(P);
    P = permute(P,[1 2 4 3]);
    sz = size(P);
    P = reshape(P,[prod(sz(1:3)) sz(4)]);
    P = P';
    [U,E] = PCA_FIT_FULL_Tws(P,7);
    pC = PCA_REPROJ_T(P,E,U);


    %%
    % init the arborist with rules for building trees
    options = statset('Display','iter');
    % creaete extract function to give to the arborist
    filter = fspecial('gaussian',[21 21],5);
    extractFunc = @(X)X;
    % create tree growth rules for the arborist
    suggestNumberOfClusterFunction = @(data,level,treePara)treePara(1)*(level<=treePara(2)) + 1*(level>treePara(2));
    % create cluster parameter generating function for arborist
    clusterFunctionFunctionGenerator = @(X,K)fitgmdist(X,K,'Options',options,'Start','plus','RegularizationValue',0.0001);
    % feature selection function
    idxSelectorFunction = @(X,L)logical(ones(1,subDims));
    % spec the tree parameters
    maxBranch = 3;
    maxDepth = 3;
    % build the arborist
    jA = arborist(suggestNumberOfClusterFunction,clusterFunctionFunctionGenerator,extractFunc,maxDepth,maxBranch,idxSelectorFunction);
    subDims = 4;
    %% plant the trees
    forest = jA.plantTrees(pC(1:subDims,:)',10);
    %%
    str = 1;
    stp = str + oneSZ - 1;
    for img = 1:30
        oneSZ = [452*452];
       
        subData = pC(:,str:stp)';
        str = stp + 1;
        stp = str + oneSZ - 1;


        [k] = forest.clusterData(subData);
        for e = 1:size(k,2)
            ki = reshape(k(:,e),[452 452]);
            kirgb = label2rgb(ki');
            imshow(kirgb,[]);
            drawnow
            waitforbuttonpress
        end
    end

%}