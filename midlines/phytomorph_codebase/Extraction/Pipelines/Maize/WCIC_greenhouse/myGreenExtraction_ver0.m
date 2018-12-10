function [GMModel] = myGreenExtraction_ver0(fileList,toFunc)
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
        tmpo = imread(fileList{e});
        tmp = toFunc(double(tmpo)/255);
        tmp = imfilter(tmp,fspecial('disk',11),'replicate');
        tmp = imresize(tmp,alpha);
        gtmp = rgb2gray(tmpo);
        bright(e) = mean(gtmp(:));
        %tmpS = decorrstretch(tmp);
        fprintf(['done with:' num2str(e) '\n'])
        I(:,:,:,e) = tmp;
    end
    rmidx =  bright < 115;
    I(:,:,:,rmidx) = [];
    I = permute(I,[1 2 4 3]);
    sz = size(I);
    I = reshape(I,[prod(sz(1:3)) sz(4)]);
    
    % fit model
    SKIP = 1;
    GMModel = fitgmdist(255*I(1:SKIP:end,:),3,'RegularizationValue',0.5,'Start','plus');
    
    
    
    %{
    % apply model
    for e = 1:numel(fileList)
        tmp = imread(fileList{e});
        tmp = toFunc(double(tmp)/255)*255;
        
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

   
    
    B = [];
    for e = 1:30
        B(:,:,:,e) = double(imread(fileList{gidx{1}(e)}))/255;
    end


    filterFunc = @(X)batchMatch(X,B);
    GMM = myGreenExtraction_ver0(fileList,filterFunc);
    [MASK] = applyGMM(fileList,filterFunc,GMM);

    cnt = 1;
    for e = 1:numel(fileList)
        if ~isempty(strfind(fileList{e},'left'))
            msk = MASK(:,:,e) == 3;
            msk = bwareaopen(msk,400);
            msk = bwlarge(msk);
            msk = imresize(msk,.25);
            msk = imopen(logical(msk),strel('disk',5,0));


            tmpG = rgb2gray(imread(fileList{e}));
            tmpGs = stdfilt(edge(tmpG),getnhood(strel('disk',21)));
            tmpGs = tmpGs.*double((MASK(:,:,e) == 1));
            [centers, radii,metric] = imfindcircles(tmpGs,[150 250],'Sensitivity',1);
            
            tmpGsM = tmpGs > graythresh(tmpGs);
            [centers, radii,metric] = imfindcircles(tmpGsM,[150 250],'Sensitivity',1);

            
            [centers, radii,metric] = imfindcircles(tmpGsM.*tmpGs,[150 250],'Sensitivity',1);

            P = sum(tmpGs,1);
            E = edge(tmpG);


            msk = imerode(msk,strel('disk',3));
            L = imfilter(double(msk),fspecial('disk',21),'replicate');
        
            POTS = imdilate(L,strel('disk',51)) == L;

            POTS = imresize(POTS,4);

            msk = logical(imresize(msk,4));
            msk = bwlarge(msk);

            %tmp = sum(MASK(:,:,e) == 1,1) > .35*size(MASK,1);
            tmp = sum(MASK(:,:,e) == 1,1);
            tmp = tmp > 100;
            K = repmat(tmp,[size(MASK,1) 1]);
            POTS = K.*imdilate(POTS,strel('disk',51,0));
            tmpI = flattenMaskOverlay(imread(fileList{e}),logical(K),.25,'b');
            tmpI = flattenMaskOverlay(tmpI,logical(POTS),.25,'g');
        
            
            
            %out(:,:,:,cnt) = flattenMaskOverlay(imread(fileList{e}),logical(msk));
            imshow(out(:,:,:,cnt),[]);
            drawnow
            dbm(cnt) = sum(msk(:));
            dn(cnt) = kidx(e);
            cnt = cnt + 1;
        end
    end


    cnt = 1;
    for e = 1:numel(fileList)
        if ~isempty(strfind(fileList{e},'left'))
            msk = MASK(:,:,e) == 3;
            msk = bwareaopen(msk,400);
            msk = bwlarge(msk);
            msk = imresize(msk,.25);
            msk = imopen(logical(msk),strel('disk',5,0));
            vec(cnt,:) = mean(msk,1);

            tmp = sum(MASK(:,:,e) == 1,1) > .75*size(msk,1);
            K = repmat(tmp,[size(MASK,1),1) 1]);
            MOV(:,:,:,cnt) = imread(fileList{e});
            cnt = cnt + 1;
        end
    end

    for e = 1:size(MOV,4)
        GR(:,:,e) = rgb2gray(MOV(:,:,:,e));
    end

    for e = 1:size(out,4)
        imwrite(out(:,:,:,e),['/mnt/spaldingdata/Edgar/WCIC_greenhouse/camera_ver0/frameStack/' num2str(e) '.tif']);
    end

    writerObj = VideoWriter('/mnt/spaldingdata/Edgar/WCIC_greenhouse/camera_ver0/frameStack/stack1.avi');
    open(writerObj);
    for e = 2:size(out,4)
        frame = im2frame(imresize(out(:,:,:,e),.35));
        writeVideo(writerObj,frame);
    end
    close(writerObj);
    sdbm = imfilter(dbm(2:end),fspecial('average',[1 11]),'replicate');
    csvwrite('/mnt/spaldingdata/Edgar/WCIC_greenhouse/camera_ver0/digital_bioMass.csv.',[sdbm'dbm(2:end)']);

    %batchMatch(fileList(gidx{3}(10)),fileList(gidx{1}(2:10)));

    for u = [1 3];
        Model{u} = myGreenExtraction(fileList(kidx==u));
    end
%}