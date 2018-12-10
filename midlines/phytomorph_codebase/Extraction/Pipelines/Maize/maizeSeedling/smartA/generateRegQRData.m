function [opTable] = generateRegQRData(opTable,GMModelFILE)

    didx = find(strcmp(opTable.type,'qr_image'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % build or load cluster
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(GMModelFILE) 
        %% collect color space
        CS = [];
        for e = 1:min(numel(midx),100)
            I = imread(filelist{e});
            sz = size(I);
            CI = reshape(I,[prod(sz(1:2)) sz(3)]);
            CS = [CS;CI];
            e
        end

        %%
        flag = true;
        while flag
            try
                %plot3(CS(1:1000:end,1),CS(1:1000:end,2),CS(1:1000:end,3),'.')
                GMModel = fitgmdist(double(CS(1:100:end,:)),3);
                flag = false;
            catch
                flag = true;
            end
        end
    else
        load('/mnt/tetra/nate/modelData/GMM.mat');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% make pipe blocks
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    W = 15;
    w = 4;
    B = zeros(2*W+1);
    B((end-1)/2+1,(end-1)/2+1) = 1;
    B((end-1)/2-w+1:(end-1)/2+w+1,1:((end-1)/2+1+w)) = 1;
    imshow(B,[]);

    BB = [];
    BB(:,:,1) = B;
    for e = 2:4
        BB(:,:,e) = imrotate(BB(:,:,e-1),90);
        imshow(BB(:,:,e),[]);
        drawnow
    end
    BL = [];
    cnt = 1;
    for e = 2:4
        v = nchoosek(1:4,e);
        for k = 1:size(v,1)
            newB = zeros(size(B));
            for i = 1:size(v,2)
                newB = newB | (BB(:,:,v(k,i)) == 1);
            end
            newB = newB .* sum(newB(:)).^-1 + e*.00001;
            BL(:,:,cnt) = newB;
            imshow(BL(:,:,cnt),[])
            drawnow
            cnt = cnt + 1;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    close all
    newM = [];
    newI = [];
    registeredQRpath = '/mnt/tetra/nate/seedlingImageParts/qrCodes/regQRcode/';



    isModelBuilt = 0;
    for f = 1:numel(didx)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the file name, frame key and read image
        fileName = opTable.data{didx(f)};
        frameKey = opTable.frameKey(didx(f));
        I = imread(fileName);
        sz = size(I);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% cluster color space
        CI = reshape(I,[prod(sz(1:2)) sz(3)]);
        T = GMModel.cluster(double(CI));
        T = reshape(T,sz(1:2));
        T = bwlarge(T==2);
        T = imclose(T,strel('square',7));
        R = [];
        for c = 1:size(BL,3)
            R(:,:,c) = imfilter(double(T),double(BL(:,:,c)),'replicate');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [IDXP,IDX] = max(R,[],3);
        eT = imerode(T,strel('square',5));
        mIDX = IDX.*eT;
        toSearchFor = [1 3 4 6 7 8 9 10];
        pt = [];
        for e = 1:numel(toSearchFor)
            ff = bwlarge(mIDX == toSearchFor(e));
            fidx = find(ff);
            [mm,midx] = max(IDXP(fidx));
            midx = find(IDXP(fidx) == mm);
            AL = [];
            [AL(:,1) AL(:,2)] = ind2sub(size(T),fidx(midx));
            pt(e,:) = mean(AL,1);
        end
        imshow(I,[]);
        hold on
        plot(pt(:,2),pt(:,1),'k*')
        drawnow
        hold on
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for e = 1:size(pt,1)
            plot(pt(e,2),pt(e,1),'k*');
            text(pt(e,2)+3,pt(e,1),num2str(e));
        end

        checkBOX = pt([8 6 4 2],:);
        qrBOX = pt([3 5 8 7],:);
        humanBOX = pt([5 1 7 6],:);

        plot(checkBOX(1:2,2),checkBOX(1:2,1),'k')
        plot(checkBOX(3:4,2),checkBOX(3:4,1),'k')
        plot(checkBOX([1 3],2),checkBOX([1 3],1),'k')
        plot(checkBOX([2 4],2),checkBOX([2 4],1),'k')

        plot(qrBOX(1:2,2),qrBOX(1:2,1),'k')
        plot(qrBOX(3:4,2),qrBOX(3:4,1),'k')
        plot(qrBOX([1 3],2),qrBOX([1 3],1),'k')
        plot(qrBOX([2 4],2),qrBOX([2 4],1),'k')

        plot(humanBOX(1:2,2),humanBOX(1:2,1),'k')
        plot(humanBOX(3:4,2),humanBOX(3:4,1),'k')
        plot(humanBOX([1 3],2),humanBOX([1 3],1),'k')
        plot(humanBOX([2 4],2),humanBOX([2 4],1),'k')
        drawnow
        pause(.1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        opTable = insertIntoMetaTable(opTable,checkBOX,'checkPolygon',frameKey);
        opTable = insertIntoMetaTable(opTable,qrBOX,'qrPolygon',frameKey);
        opTable = insertIntoMetaTable(opTable,humanBOX,'humanPolygon',frameKey);
        opTable = insertIntoMetaTable(opTable,pt,'framePoints',frameKey);


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        pt = fliplr(pt);

        close all
        if f ~= 1 | isModelBuilt
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % transform to the model - the first image
            tform = fitgeotrans(pt,fixedPoints,'similarity');
            tmpI = imwarp(I,tform,'OutputView',outP);
            newI(:,:,:,f) = tmpI;
            newM(:,:,f) = imwarp(T,tform,'OutputView',outP);
            


            toV = mean(newI/255,4);
            toM = mean(newM,3);
            out = flattenMaskOverlay(toV,toM > .3,.8,'g');
            imshow(out,[]);
            drawnow
            hold off
            pause(.4)
            %waitforbuttonpress
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % init the model to the first image in the list
            % set the fixed points to the first image
            fixedPoints = pt;
            outSize = size(I);
            outP = imref2d(outSize(1:2));
            newI(:,:,:,1) = I;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end

        %{
        if isModelBuilt
            imwrite(newI(:,:,:,i),[registeredQRpath nm '.tif']);
        end
        %}
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    here = 1;
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% build model
    toV = mean(newI/255,4);
    toM = mean(newM,3);
    frameModel = bsxfun(@times,(toM > .5),toV);
    T = (toM > .5);
    R = [];
    for c = 1:size(BL,3)
        R(:,:,c) = imfilter(double(T),double(BL(:,:,c)),'replicate');
    end
    [IDXP,IDX] = max(R,[],3);
    eT = imerode(T,strel('square',5));
    mIDX = IDX.*eT;
    %imshow(mIDX,[]);
    toSearchFor = [1 3 4 6 7 8 9 10];
    pt = [];
    for e = 1:numel(toSearchFor)
        ff = bwlarge(mIDX == toSearchFor(e));
        fidx = find(ff);
        [mm,midx] = max(IDXP(fidx));
        midx = find(IDXP(fidx) == mm);
        AL = [];
        [AL(:,1) AL(:,2)] = ind2sub(size(T),fidx(midx));
        pt(e,:) = mean(AL,1);
    end

    close all
    imshow(frameModel,[]);
    hold on
    for e = 1:size(pt,1)
        plot(pt(e,2),pt(e,1),'k*');
        text(pt(e,2)+3,pt(e,1),num2str(e));
    end


    checkBOX = pt([8 6 4 2],:);
    qrBOX = pt([3 5 8 7],:);
    humanBOX = pt([5 1 7 6],:);



    [a(1)] = measurePointBox(checkBOX);
    [a(2)] = measurePointBox(qrBOX);
    [a(3)] = measurePointBox(humanBOX);


    toV = imrotate(toV,mean(a)*180/pi);
    toM = imrotate(toM,mean(a)*180/pi);



    frameModel = bsxfun(@times,(toM > .5),toV);
    imshow(toV,[]);
    drawnow
    T = (toM > .5);
    R = [];
    for c = 1:size(BL,3)
        R(:,:,c) = imfilter(double(T),double(BL(:,:,c)),'replicate');
    end
    [IDXP,IDX] = max(R,[],3);
    eT = imerode(T,strel('square',5));
    mIDX = IDX.*eT;
    %imshow(mIDX,[]);
    toSearchFor = [1 3 4 6 7 8 9 10];
    pt = [];
    for e = 1:numel(toSearchFor)
        ff = bwlarge(mIDX == toSearchFor(e));
        fidx = find(ff);
        [mm,midx] = max(IDXP(fidx));
        midx = find(IDXP(fidx) == mm);
        AL = [];
        [AL(:,1) AL(:,2)] = ind2sub(size(T),fidx(midx));
        pt(e,:) = mean(AL,1);
    end
    a = [];
    close all
    imshow(frameModel,[]);
    hold on
    for e = 1:size(pt,1)
        plot(pt(e,2),pt(e,1),'k*');
        text(pt(e,2)+3,pt(e,1),num2str(e));
    end
    checkBOX = pt([8 6 4 2],:);
    qrBOX = pt([3 5 8 7],:);
    humanBOX = pt([5 1 7 6],:);
    [a(1)] = measurePointBox(checkBOX);
    [a(2)] = measurePointBox(qrBOX);
    [a(3)] = measurePointBox(humanBOX);
    fixedPoints = fliplr(pt);

    outSize = size(frameModel);
    outP = imref2d(outSize(1:2));

    %}
end


