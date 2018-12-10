function [GMModel] = myGreenExtraction_ver0(fileList)
%{    
    tmp = imread(fileList{e});
    tmp = imfilter(tmp,fspecial('disk',11),'replicate');
    tU = [10 230 10];
    tmpS = decorrstretch(tmp,'TargetMean',tU);
%}

    alpha = .30;
    N = 40;
    I = [];
    for e = 1:min(N,numel(fileList))
        tmp = imread(fileList{e});
        tmp = imfilter(tmp,fspecial('disk',11),'replicate');
        tmp = imresize(tmp,alpha);
        gtmp = rgb2gray(tmp);
        bright(e) = mean(gtmp(:));
        
       
        %tmpS = decorrstretch(tmp);
        
        I(:,:,:,e) = tmp;
    end
    rmidx =  bright < 50;
    I(:,:,:,rmidx) = [];
    I = permute(I,[1 2 4 3]);
    sz = size(I);
    I = reshape(I,[prod(sz(1:3)) sz(4)]);
    
    % fit model
    SKIP = 100;
    GMModel = fitgmdist(I(1:SKIP:end,:),4,'RegularizationValue',0.6,'Start','plus');
    
    
    
    %{
    % apply model
    for e = 1:numel(fileList)
        tmp = imread(fileList{e});
        
        tmp = imfilter(tmp,fspecial('disk',11),'replicate');
        %tmp = decorrstretch(tmp);
        
        sz = size(tmp);
        tmpo = tmp;
        tmp = double(reshape(tmp,[size(tmp,1)*size(tmp,2) size(tmp,3)]));
        idx = cluster(GMModel,tmp);
        idx = reshape(idx,sz(1:2));
        MASK = idx == 1;
        RGB = label2rgb(idx);
        imshow([RGB tmpo],[]);
        drawnow
    end
    %}
    
    done = 1;
end

%{
    FilePath = '/home/nate/Downloads/WCICgreenhouse/'
    %FilePath = '/home/nate/Downloads/right/';
    FileList = {};
    FileExt = {'tiff','TIF','JPG'};
    verbose = 1;
    SET = sdig(FilePath,FileList,FileExt,verbose);
    [fileList] = stackFileList(SET);
    H = getBulkHistogram(fileList);
    [hS hC hU hE hL hERR hLAM] = PCA_FIT_FULL(H,3);
    kidx = kmeans(hC,3);
    UQ = unique(kidx);
    for u = 1:numel(UQ)
        fidx = find(kidx==UQ(u));
        for k = 1:min(numel(fidx),10)
            I = imread(fileList{fidx(k)});
            imshow(I,[]);
            title(num2str(u))
            drawnow
        end
    end
    
    for u = 1:3
        gidx{u} = find(kidx==u);
    end

    [MASK] = applyGMM(fileList,GMModel)

    

    myGreenExtraction(fileList,



    %batchMatch(fileList(gidx{3}(10)),fileList(gidx{1}(2:10)));

    for u = [1 3];
        Model{u} = myGreenExtraction(fileList(kidx==u));
    end
%}