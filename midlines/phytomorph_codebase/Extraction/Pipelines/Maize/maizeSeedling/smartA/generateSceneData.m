function [] = generateSceneData(fileName,BackgroundGMM,nonBackgroundGMM,plantBOXresize,conetainerBOXresize,qrBOXresize,plantBufferWidth,bioDataBuffer,coneTainerModelDepth,oPath)
    % bioDataBuffer : the extra amount to clip off from the bottom of the QR code
   
    
    objectTable = table;
    baseLocation = '/mnt/tetra/nate/seedlingImageParts/';
    
    
    % get the file parts
    [~,nm] = fileparts(fileName);
    % read the image
    I = double(imread(fileName))/255;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % handle image split
    % cut the image
    [dataStrip,bioStrip,cropLine,msg,qrCropBox] = splitMaizeSeedlingImage(I,bioDataBuffer);
    % buffer up the QR code to include blue background
    qrCropBox(1) = qrCropBox(1) - 50;
    qrCropBox(2) = qrCropBox(2) - 50;
    qrCropBox(3) = qrCropBox(3) + 100;
    qrCropBox(4) = qrCropBox(4) + 100;
    
    
    
    RAW_sceneStruct.bioDataDivideLocation = size(bioStrip,1);
    RESIZED_sceneStruct.bioDataDivideLocation = size(bioStrip,1);
    RAW_sceneStruct.bioDataBuffer = bioDataBuffer;
    RESIZED_sceneStruct.bioDataBuffer = bioDataBuffer;
    
    
    
    % create masks for background and conetainers
    [M,LABEL] = createMaskSets(bioStrip,BackgroundGMM,nonBackgroundGMM);
    % clean the conetainers
    [coneTainer,coneTainerCropBox] = processMaskSet(I,M,100);
    RAW_sceneStruct.coneTainerMask = coneTainer;
    RAW_sceneStruct.backgroundMask = M(:,:,1);
    RAW_sceneStruct.nonBackgroundMask = M(:,:,3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % handle QR data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % gather the qr code block for modeling
    tmpQR = imcrop(I,qrCropBox);
    qrMASK = zeros(size(dataStrip,1),size(dataStrip,2));
    qrMASK(round(qrCropBox(2):qrCropBox(2)+qrCropBox(4)),round(qrCropBox(1):qrCropBox(1)+qrCropBox(3))) = 1;
    RAW_sceneStruct.qrMask = qrMASK;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% WRITE QR RAW
    RAW_sceneStruct.qrObject.imageData = tmpQR;
    RAW_sceneStruct.qrObject.cropBox = qrCropBox;
    %%%%%% WRITE QR RESIZED
    RESIZED_sceneStruct.qrObject.imageData = imresize(tmpQR,qrBOXresize);
    RESIZED_sceneStruct.qrObject.cropBox = qrCropBox;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% WRITE QR to object list and table
    qrFileName = [baseLocation 'qrCodes' filesep 'qrRaw' filesep nm '.tif'];
    pt = size(objectTable,1);
    objectTable{pt+1,'type'} = {'qrObjectRaw'};
    objectTable{pt+1,'imageLocation'} = {qrFileName};
    for co = 1:4
        objectTable{pt+1,['boundingBox' num2str(co)]} = qrCropBox(co);
    end
    imwrite(tmpQR,qrFileName);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % handle background modeling
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    wholeBK = [~qrMASK;M(:,:,1)];
    looseBK = M(:,:,1);
    looseBK = imerode(looseBK,strel('disk',5,0));
    TOT = sum(looseBK,1);
    mI = bsxfun(@times,bioStrip,looseBK);
    BK1 = [];
    for k = 1:size(mI,3)
        BK1(k,:) = sum(mI(:,:,k),1).*TOT.^-1;
        BK1(k,:) = fillmissing(BK1(k,:),'nearest');
        BK1(k,:) = imfilter(BK1(k,:),fspecial('average',[1 401]),'replicate');
    end
    SIM1 = shiftdim(BK1',-1);
    SIM1 = repmat(SIM1,[size(bioStrip,1) 1 1]);
    
    
    DELTABK = bsxfun(@times,mI - SIM1,looseBK);
    DELTABK(find(isnan(DELTABK))) = 0;
    TOT = sum(looseBK,2);
    BK2 = [];
    for k = 1:size(mI,3)
        BK2(k,:) = squeeze(sum(DELTABK(:,:,k),2).*TOT.^-1);
        BK2(k,:) = imfilter(BK2(k,:),fspecial('average',[1 401]),'replicate');
    end
    
    for k = 1:size(BK1,1)
         iBK1(k,:) = interp1(BK1(k,:),linspace(1,size(BK1,2),5000));
         iBK2(k,:) = interp1(BK2(k,:),linspace(1,size(BK2,2),2500));
    end
    RESIZED_sceneStruct.iBK1 = iBK1;
    RESIZED_sceneStruct.iBK2 = iBK2;

    
    %%%%%% WRITE background to object list and table
    tmpFileName = [baseLocation 'backgrounds' filesep 'h' filesep nm '.tif'];
    pt = size(objectTable,1);
    objectTable{pt+1,'type'} = {'background_h'};
    objectTable{pt+1,'imageLocation'} = {tmpFileName};
    tmpData = shiftdim(iBK1',-1);
    imwrite(tmpData,tmpFileName);
    
    
    %%%%%% WRITE background to object list and table
    tmpFileName = [baseLocation 'backgrounds' filesep 'v' filesep nm '.tif'];
    pt = size(objectTable,1);
    objectTable{pt+1,'type'} = {'background_v'};
    objectTable{pt+1,'imageLocation'} = {tmpFileName};
    tmpData = shiftdim(iBK2',-1);
    tmpData = permute(tmpData,[2 1 3]);
    imwrite(tmpData,tmpFileName);
    
    
    
    
    % simlulate background
    simBackground1 = repmat(imresize(shiftdim(iBK1',[-1]),[1 size(I,2)]),[size(I,1) 1 1]);
    simBackground2 = repmat(imresize(permute(shiftdim(iBK2',[-1]),[2 1 3]),[size(I,1) 1]),[1 size(I,2) 1]);
    simBackground = simBackground1 + simBackground2;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % locate other objects in top
    deltaI = simBackground - I;
    deltaI = sum(deltaI.^2,3).^.5;
    deltaI(qrCropBox(2):qrCropBox(2)+qrCropBox(4)-1,qrCropBox(1):qrCropBox(1)+qrCropBox(3)-1,:) = 0;
    deltaIData = deltaI(1:cropLine,:,:);
    extraObj = deltaIData > graythresh(deltaIData);
    extraObj = bwareaopen(extraObj,200);
    extraObj = imfill(extraObj,'holes');
    extraObj = imclearborder(extraObj);
    extraObj = bwareaopen(extraObj,5000);
    eR = regionprops(logical(extraObj),'Area','BoundingBox');
    if ~isempty(eR)
        for ex = 1:numel(eR)
            eR(ex).BoundingBox(1) = eR(ex).BoundingBox(1) - 50;
            eR(ex).BoundingBox(2) = eR(ex).BoundingBox(2) - 50;
            eR(ex).BoundingBox(3) = eR(ex).BoundingBox(3) + 50;
            eR(ex).BoundingBox(4) = eR(ex).BoundingBox(4) - 50;
            extraObject = imcrop(dataStrip,eR(ex).BoundingBox);
            
            %%%%%% WRITE QR to object list and table
            tmpFileName = [baseLocation 'extraDataObjects' filesep nm '_' num2str(ex) '.tif'];
            pt = size(objectTable,1);
            objectTable{pt+1,'type'} = {'extraDataObject'};
            objectTable{pt+1,'imageLocation'} = {tmpFileName};
            for co = 1:4
                if co == 2
                    offset = cropLine;
                else 
                    offset = 0;
                end
                objectTable{pt+1,['boundingBox' num2str(co)]} = qrCropBox(co) + offset;
            end
            imwrite(extraObject,tmpFileName);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate QR SIM
    QRSIM = zeros(size(I));
    QRSIM = simBackground;
    qrLoc = [1000 100];
    QR_blend = zeros(size(I));
    QR_blend(qrLoc(2):qrLoc(2)+size(tmpQR,1)-1,qrLoc(1):qrLoc(1)+size(tmpQR,2)-1,:) = 1;
    QR_blend = imerode(QR_blend,strel('disk',21,0));
    QR_blend = imfilter(QR_blend,fspecial('gaussian',[71 71],21),'replicate');
    QRSIM(qrLoc(2):qrLoc(2)+size(tmpQR,1)-1,qrLoc(1):qrLoc(1)+size(tmpQR,2)-1,:) = tmpQR;
    blendSIM = QR_blend.*QRSIM + (1-QR_blend).*simBackground;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    
    plantBox = [];
    % create mod crop boxes for each container
    for r = 1:numel(coneTainerCropBox)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create mod crop box for stable modeling of fixed height
        coneTainerCropBoxMOD(r).BoundingBox = coneTainerCropBox(r).BoundingBox;
        coneTainerCropBoxMOD(r).BoundingBox(4) = coneTainerModelDepth;
        
        coneTainerCropBoxMODLONG(r).BoundingBox = coneTainerCropBox(r).BoundingBox;
        coneTainerCropBoxMODLONG(r).BoundingBox(4) = size(bioStrip,1) - coneTainerCropBoxMODLONG(r).BoundingBox(2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % crop the whole height conetainer
        coneTemp = imcrop(bioStrip,coneTainerCropBoxMODLONG(r).BoundingBox);
        coneLAB = imcrop(LABEL,coneTainerCropBoxMODLONG(r).BoundingBox);
        coneMASK = coneLAB == 3 | coneLAB == 1 | coneLAB == 4;
        %%%%%% WRITE container to object list and table
        tmpFileName = [baseLocation 'conetainers' filesep 'whole' filesep nm '_' num2str(r) '.tif'];
        pt = size(objectTable,1);
        objectTable{pt+1,'type'} = {'conetainer_whole'};
        objectTable{pt+1,'imageLocation'} = {tmpFileName};
        for co = 1:4
            
            if co == 2
                offset = cropLine;
            else 
                offset = 0;
            end
            objectTable{pt+1,['boundingBox' num2str(co)]} = coneTainerCropBoxMODLONG(r).BoundingBox(co) + offset;
            
        end
        imwrite(coneTemp,tmpFileName);
        
        %%%%%% WRITE container to object list and table
        tmpFileName = [baseLocation 'conetainers' filesep 'whole_masked' filesep nm '_' num2str(r) '.tif'];
        pt = size(objectTable,1);
        objectTable{pt+1,'type'} = {'conetainer_whole_masked'};
        objectTable{pt+1,'imageLocation'} = {tmpFileName};
        for co = 1:4
            
            if co == 2
                offset = cropLine;
            else 
                offset = 0;
            end
            objectTable{pt+1,['boundingBox' num2str(co)]} = coneTainerCropBoxMODLONG(r).BoundingBox(co) + offset;
            
        end
        coneTempMASKED = bsxfun(@times,coneTemp,coneMASK);
        imwrite(coneTempMASKED,tmpFileName);

        
        
        
        
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % crop the fixed height conetainer
        coneTemp = imcrop(bioStrip,coneTainerCropBoxMOD(r).BoundingBox);
        % store conetainer images for stable
        RAW_sceneStruct.conetainerBox(r).IMAGE = coneTemp;
        RESIZED_sceneStruct.conetainerBox(r).IMAGE = imresize(coneTemp,conetainerBOXresize);
        
        %%%%%% WRITE container to object list and table
        tmpFileName = [baseLocation 'conetainers' filesep 'clipped' filesep nm '_' num2str(r) '.tif'];
        pt = size(objectTable,1);
        objectTable{pt+1,'type'} = {'conetainer_clipped'};
        objectTable{pt+1,'imageLocation'} = {tmpFileName};
        for co = 1:4
            
            
            if co == 2
                offset = cropLine;
            else 
                offset = 0;
            end
            objectTable{pt+1,['boundingBox' num2str(co)]} = coneTainerCropBoxMOD(r).BoundingBox(co) + offset;
            
            
        end
        imwrite(coneTemp,tmpFileName);
        
        
        plantBox(r).BoundingBox = coneTainerCropBox(r).BoundingBox;
        plantBox(r).BoundingBox(2) = 1;
        plantBox(r).BoundingBox(4) = coneTainerCropBox(r).BoundingBox(2)-1;
        plantBox(r).BoundingBox(1) = plantBox(r).BoundingBox(1) - plantBufferWidth;
        plantBox(r).BoundingBox(3) = plantBox(r).BoundingBox(3) + 2*plantBufferWidth;

        
        
        %%%%%%%%%%%%%%%%%%
        % store centroid of conetainer
        RAW_sceneStruct.conetainerBox(r).Centroid = coneTainerCropBox(r).Centroid;
        RESIZED_sceneStruct.conetainerBox(r).Centroid = coneTainerCropBox(r).Centroid;
        %%%%%%%%%%%%%%%%%%
        % store index number of plant and container
        RAW_sceneStruct.plantBox(r).coneTainerLocation = r;
        RESIZED_sceneStruct.plantBox(r).coneTainerLocation = r;
        RAW_sceneStruct.conetainerBox(r).coneTainerLocation = r;
        RESIZED_sceneStruct.conetainerBox(r).coneTainerLocation = r;
        %%%%%%%%%%%%%%%%%%
        % store bounding boxes for plants
        RAW_sceneStruct.plantBox(r).BoundingBox = plantBox(r).BoundingBox;
        RESIZED_sceneStruct.plantBox(r).BoundingBox = plantBox(r).BoundingBox;
        %%%%%%%%%%%%%%%%%%
        % store bounding boxes for conetainers
        RAW_sceneStruct.conetainerBox(r).BoundingBox = coneTainerCropBox(r).BoundingBox;
        RESIZED_sceneStruct.conetainerBox(r).BoundingBox = coneTainerCropBox(r).BoundingBox;
        %%%%%%%%%%%%%%%%%%
        % store bounding boxes for conetainers
        RAW_sceneStruct.conetainerBox(r).MODBoundingBox = coneTainerCropBoxMOD(r).BoundingBox;
        RESIZED_sceneStruct.conetainerBox(r).MODBoundingBox = coneTainerCropBoxMOD(r).BoundingBox;
        %%%%%%%%%%%%%%%%%%
        % store masks
        w = imcrop(M(:,:,3),plantBox(r).BoundingBox);
        m = bwlarge(w);
        RAW_sceneStruct.plantBox(r).MASK = m;
        RESIZED_sceneStruct.plantBox(r).MASK = imresize(m,plantBOXresize);
        
        
        %%%%%% WRITE container to object list and table
        tmpFileName = [baseLocation 'plants' filesep 'masks' filesep nm '_' num2str(r) '.tif'];
        pt = size(objectTable,1);
        objectTable{pt+1,'type'} = {'plant_mask'};
        objectTable{pt+1,'imageLocation'} = {tmpFileName};
        for co = 1:4
            
            
            if co == 2
                offset = cropLine;
            else 
                offset = 0;
            end
            objectTable{pt+1,['boundingBox' num2str(co)]} = plantBox(r).BoundingBox(co) + offset;
            
        end
        imwrite(m,tmpFileName);
        
        
        %%%%%%%%%%%%%%%%%%
        % store images
        w = imcrop(bioStrip,plantBox(r).BoundingBox);
        RAW_sceneStruct.plantBox(r).IMAGE = w;
        RESIZED_sceneStruct.plantBox(r).IMAGE = imresize(w,plantBOXresize);
        
        
        %%%%%% WRITE container to object list and table
        tmpFileName = [baseLocation 'plants' filesep 'color' filesep nm '_' num2str(r) '.tif'];
        pt = size(objectTable,1);
        objectTable{pt+1,'type'} = {'plant_color'};
        objectTable{pt+1,'imageLocation'} = {tmpFileName};
        for co = 1:4
            if co == 2
                offset = cropLine;
            else 
                offset = 0;
            end
            objectTable{pt+1,['boundingBox' num2str(co)]} = plantBox(r).BoundingBox(co) + offset;
        end
        imwrite(w,tmpFileName);
        
        
        
        %%%%%% WRITE container to object list and table
        w = bsxfun(@times,w,m);
        tmpFileName = [baseLocation 'plants' filesep 'colorMasked' filesep nm '_' num2str(r) '.tif'];
        pt = size(objectTable,1);
        objectTable{pt+1,'type'} = {'plant_color_masked'};
        objectTable{pt+1,'imageLocation'} = {tmpFileName};
        for co = 1:4
            if co == 2
                offset = cropLine;
            else 
                offset = 0;
            end
            objectTable{pt+1,['boundingBox' num2str(co)]} = plantBox(r).BoundingBox(co) + offset;
        end
        imwrite(w,tmpFileName);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    
    

    %generateLabledScene(I,RAW_sceneStruct)
    tableName = [baseLocation 'dataTables' filesep nm '.txt']; 
    writetable(objectTable,tableName);
    
    [~,nm] = fileparts(fileName);
    if ~isempty(oPath)
        save([oPath nm '.mat'],'RESIZED_sceneStruct','RAW_sceneStruct');
    end



end