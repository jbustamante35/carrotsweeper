function [opTable] = generateSceneData_ver2(opTable,BackgroundGMM,nonBackgroundGMM,plantBufferWidth,bioDataBuffer,coneTainerModelDepth,baseLocation)
    

    % bioDataBuffer : the extra amount to clip off from the bottom of the QR code
    
    
    while opTable.hasdata
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set up the run - get file name, frame and read image
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % get the frame key for the file
        [fileName,info] = opTable.read();
        % mod the start time
        opTable.modExecuteContext('state','loading',info.executeKey);
        % get the fth file name
        frameKey = info.frameKey;
        contextKey = info.executeKey;
        % get the file parts
        [~,nm] = fileparts(fileName);
        % read the image
        I = double(imread(fileName))/255;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % handle image split
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [dataStrip,bioStrip,cropLine,msg,qrCropBox] = splitMaizeSeedlingImage(I,bioDataBuffer);
        % buffer up the QR code to include blue background
        qrCropBox(1) = qrCropBox(1) - 50;
        qrCropBox(2) = qrCropBox(2) - 50;
        qrCropBox(3) = qrCropBox(3) + 100;
        qrCropBox(4) = qrCropBox(4) + 100;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % insert qr crop box
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        qrKey = opTable.put(qrCropBox,'qr_CropBox',frameKey);
        bufferSegment = [cropLine 1 cropLine size(I,2)];
        bufferKey = opTable.put(bufferSegment,'buffer_LineSegment',frameKey);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create masks for background and conetainers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [M,LABEL] = createMaskSets(bioStrip,BackgroundGMM,nonBackgroundGMM);
        % clean the conetainers
        [coneTainer,coneTainerCropBox] = processMaskSet(I,M,100);
        RAW_sceneStruct.coneTainerMask = coneTainer;
        RAW_sceneStruct.backgroundMask = M(:,:,1);
        RAW_sceneStruct.nonBackgroundMask = M(:,:,3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % handle QR data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % gather the qr code block for modeling
        tmpQR = imcrop(I,qrCropBox);
        qrMASK = zeros(size(dataStrip,1),size(dataStrip,2));
        qrMASK(round(qrCropBox(2):qrCropBox(2)+qrCropBox(4)),round(qrCropBox(1):qrCropBox(1)+qrCropBox(3))) = 1;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% WRITE QR to object list and table
        
        opTable.imwrite(tmpQR,'qrCodes',contextKey);
        opTable.imwrite(tmpQR);
        
        
        qrFileName = [baseLocation 'qrCodes' filesep 'qrRaw' filesep nm '.tif'];
        imwrite(tmpQR,qrFileName);
        qrFrame = eye(3);
        qrFrame(:,3) = [fliplr(qrCropBox(1:2)) 1]';
        [opTable,~,qrFrameKey] = insertIntoMetaTable(opTable,qrFrame,'qr_frame',frameKey);
        opTable = insertIntoMetaTable(opTable,qrFileName,'qr_image',qrFrameKey);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %{
        %insertIntoMetaTable(opTable,qrFileName,'qrImage',
        pt = size(objectTable,1);
        objectTable{pt+1,'type'} = {'qrObjectRaw'};
        objectTable{pt+1,'imageLocation'} = {qrFileName};
        for co = 1:4
            objectTable{pt+1,['boundingBox' num2str(co)]} = qrCropBox(co);
        end
        imwrite(tmpQR,qrFileName);
        %}
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
        

        %%%%%% WRITE background to object list and table
        tmpFileName = [baseLocation 'backgrounds' filesep 'h' filesep nm '.tif'];
        tmpData = shiftdim(iBK1',-1);
        imwrite(tmpData,tmpFileName);
        opTable = insertIntoMetaTable(opTable,tmpFileName,'h_background_image',frameKey);

        
        %%%%%% WRITE QR to object list and table
        tmpFileName = [baseLocation 'backgrounds' filesep 'v' filesep nm '.tif'];
        tmpData = shiftdim(iBK2',-1);
        tmpData = permute(tmpData,[2 1 3]);
        imwrite(tmpData,tmpFileName);
        opTable = insertIntoMetaTable(opTable,tmpFileName,'v_background_image',frameKey);
        
        
        % simlulate background
        simBackground1 = repmat(imresize(shiftdim(iBK1',[-1]),[1 size(I,2)]),[size(I,1) 1 1]);
        simBackground2 = repmat(imresize(permute(shiftdim(iBK2',[-1]),[2 1 3]),[size(I,1) 1]),[1 size(I,2) 1]);
        simBackground = simBackground1 + simBackground2;



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
                imwrite(extraObject,tmpFileName);
                tmpFrame = eye(3);
                tmpFrame(:,3) = [fliplr(eR(ex).BoundingBox(1:2)) 1]';
                [opTable,~,tmpFrameKey] = insertIntoMetaTable(opTable,tmpFrame,'extra_frame',frameKey);
                opTable = insertIntoMetaTable(opTable,tmpFileName,'extraDataObject_image',tmpFrameKey);
                
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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
            imwrite(coneTemp,tmpFileName);
            tmpFrame = eye(3);
            tmpFrame(:,3) = [fliplr(coneTainerCropBoxMODLONG(r).BoundingBox(1:2)) 1]';
            [opTable,~,tmpFrameKey] = insertIntoMetaTable(opTable,tmpFrame,'conetainer_frame',frameKey);
            opTable = insertIntoMetaTable(opTable,tmpFileName,'conetainer_whole_image',tmpFrameKey);
            
            
            
            
            %%%%%% WRITE container to object list and table
            tmpFileName = [baseLocation 'conetainers' filesep 'whole_masked' filesep nm '_' num2str(r) '.tif'];
            coneTempMASKED = bsxfun(@times,coneTemp,coneMASK);
            imwrite(coneTempMASKED,tmpFileName);
            tmpFrame = eye(3);
            tmpFrame(:,3) = [fliplr(coneTainerCropBoxMODLONG(r).BoundingBox(1:2)) 1]';
            [opTable,~,tmpFrameKey] = insertIntoMetaTable(opTable,tmpFrame,'conetainerMasked_frame',frameKey);
            opTable = insertIntoMetaTable(opTable,tmpFileName,'conetainer_wholemasked_image',tmpFrameKey);
           



            %{
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % crop the fixed height conetainer
            coneTemp = imcrop(bioStrip,coneTainerCropBoxMOD(r).BoundingBox);
            

            %%%%%% WRITE container to object list and table
            tmpFileName = [baseLocation 'conetainers' filesep 'clipped' filesep nm '_' num2str(r) '.tif'];
            imwrite(coneTemp,tmpFileName);
            opTable = insertIntoMetaTable(opTable,tmpFileName,'conetainer_clipped_image',frameKey);
            %}

            plantBox(r).BoundingBox = coneTainerCropBox(r).BoundingBox;
            plantBox(r).BoundingBox(2) = 1;
            plantBox(r).BoundingBox(4) = coneTainerCropBox(r).BoundingBox(2)-1;
            plantBox(r).BoundingBox(1) = plantBox(r).BoundingBox(1) - plantBufferWidth;
            plantBox(r).BoundingBox(3) = plantBox(r).BoundingBox(3) + 2*plantBufferWidth;


            %{
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
            %}
            
            
            %%%%%%%%%%%%%%%%%%
            % store masks
            w = imcrop(M(:,:,3),plantBox(r).BoundingBox);
            m = bwlarge(w);

            %%%%%% WRITE container to object list and table
            tmpFileName = [baseLocation 'plants' filesep 'masks' filesep nm '_' num2str(r) '.tif'];
            imwrite(m,tmpFileName);
            tmpFrame = eye(3);
            tmpFrame(:,3) = [fliplr(plantBox(r).BoundingBox(1:2)) 1]';
            [opTable,~,tmpFrameKey] = insertIntoMetaTable(opTable,tmpFrame,'plant_mask_frame',frameKey);
            opTable = insertIntoMetaTable(opTable,tmpFileName,'plant_mask_image',tmpFrameKey);
            
            %%%%%%%%%%%%%%%%%%
            % store images
            w = imcrop(bioStrip,plantBox(r).BoundingBox);
            
            %%%%%% WRITE container to object list and table
            tmpFileName = [baseLocation 'plants' filesep 'color' filesep nm '_' num2str(r) '.tif'];
            imwrite(w,tmpFileName);
            tmpFrame = eye(3);
            tmpFrame(:,3) = [fliplr(plantBox(r).BoundingBox(1:2)) 1]';
            [opTable,~,tmpFrameKey] = insertIntoMetaTable(opTable,tmpFrame,'plant_color_frame',frameKey);
            opTable = insertIntoMetaTable(opTable,tmpFileName,'plant_color_image',tmpFrameKey);

            %%%%%% WRITE container to object list and table
            tmpFileName = [baseLocation 'plants' filesep 'colorMasked' filesep nm '_' num2str(r) '.tif'];
            w = bsxfun(@times,w,m);
            imwrite(w,tmpFileName);
            tmpFrame = eye(3);
            tmpFrame(:,3) = [fliplr(plantBox(r).BoundingBox(1:2)) 1]';
            [opTable,~,tmpFrameKey] = insertIntoMetaTable(opTable,tmpFrame,'plant_color_masked_frame',frameKey);
            opTable = insertIntoMetaTable(opTable,tmpFileName,'plant_color_masked_image',tmpFrameKey);
           
        end




        %generateLabledScene(I,RAW_sceneStruct)
        tableName = [baseLocation 'dataTables' filesep nm '.mat']; 
        save(tableName,'opTable');
        %writetable(opTable,tableName);
        
        %{
        [~,nm] = fileparts(fileName);
        if ~isempty(oPath)
            save([oPath nm '.mat'],'RESIZED_sceneStruct','RAW_sceneStruct');
        end
        %}
    end



end