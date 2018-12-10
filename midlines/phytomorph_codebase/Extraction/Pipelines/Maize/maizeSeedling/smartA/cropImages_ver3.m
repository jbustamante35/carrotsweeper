function [oI,subI,boundingBoxes,connectedMASK,MASK,SKELETON,LEN] = cropImages_ver3(I,TOP_THRESH,cluster_Level0,cluster_Level1,cluster_Level2,cluster_Level3,kE,kU)


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read image if file name passed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ischar(I)
        tic;
        fprintf(['start::reading file \n']);
        I = double(imread(I));
        [I angle] = rectifyImage(I/255);
        fprintf(['end::reading file::' num2str(toc) '\n']);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read image if file name passed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    BUMP = 50;
    [msg,qrCropBox] = getQRcode(I);
    TOPloc = qrCropBox(2) + qrCropBox(4);
    TOPloc = TOPloc + BUMP;
    
    
    oI = I;
    
    LEN = NaN*zeros(1,6);
    metaHeader = I(1:TOPloc,:,:);
    hsvHeader = rgb2hsv(metaHeader);
    standardMask = hsvHeader(:,:,3) < .5;
    standardMask = imfill(standardMask,'holes');
    standardMask = bwlarge(standardMask);
    R = regionprops(standardMask,'BoundingBox','Area');
    ccm = [];
    if R.Area > 500000
        imaData = imcrop(metaHeader,R(1).BoundingBox);
        %imaData = imresize(imaData,3);
        imaData = uint8(255*imaData);
        imaDataRaster = rgb2lin(imaData);
        %imaDataRaster = rgb2lin(II);
        chart = esfrChart(imaDataRaster);
        [colorTable,ccm] = measureColor(chart);
        pts = chart.RegistrationPoints;
        pr = nchoosek(1:4,2);
        for e = 1:size(pr,1)
            toP = [pts(pr(e,:),:)];
            LEN(e) = norm(diff(toP,1,1));
        end
    end
    
    
    
    
    
    I(1:TOPloc,:,:) = [];
    
    
    tmpHSV = rgb2hsv(I);
    toLook = 2;
    %uS2(:,e) = mean(tmpHSV(:,:,toLook),2);
    %sS2(:,e) = std(tmpHSV(:,:,toLook),1,2);
    uS1 = mean(tmpHSV(:,:,toLook),1)';
   
    sS1 = std(tmpHSV(:,:,toLook),1,1)';
    sS1 = bindVec(sS1);
    uS1 = bindVec(uS1);
    sig = sS1.*uS1;
    ssig = imfilter(sig,fspecial('average',[501 1]),'circular');
    sv = sort(ssig);
    
   
    ssig(isnan(ssig)) = min(ssig(:));
    BV = zeros(1,size(I,2));
   
    
    ssig = bindVec(ssig);
    
    
    G = rgb2gray(I);
    uG = mean(imcomplement(G),1);
    uG = imfilter(uG,fspecial('average',[1 501]),'circular');
    uG = bindVec(uG);
    ssig = (ssig' + uG)';

    ssig = bindVec(ssig);
    ssig = imfilter(ssig',fspecial('average',[1 501]),'replicate');

    pk = imdilate(ssig,strel('disk',501)) == ssig;
    val = imerode(ssig,strel('disk',501)) == ssig;

    peaksPK = find(pk);
    peaksVAL = find(val);

    peaksVAL(peaksVAL < peaksPK(1)) = [];
    peaksVAL(peaksVAL > peaksPK(end)) = [];

    threshV = graythresh(ssig);
    BV = (ssig < threshV);
    BV = imclose(BV,strel('disk',40,0));
    BV = logical([zeros(size(BV));BV;zeros(size(BV))]);
    BV = imclearborder(BV);
    BV = bwlarge(BV,2);
    
    

    BV = zeros(1,size(BV,2));
    BV(peaksVAL) = 1;

    R = regionprops(logical(BV),'Centroid');


    %BV = BV(2,:);
    BV = zeros(size(BV));
    for e = 1:numel(R)
        BV(1,round(R(e).Centroid(1))) = 1;
    end
    BV(1) = 1;
    BV(end) = 1;
    
    pts = find(BV);
    
    for b = 1:(numel(pts)-1)
        CB{b} = [pts(b) TOPloc pts(b+1)-pts(b) size(oI,1)-TOPloc];
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for b = 1:numel(CB)
        J = imcrop(oI,CB{b});
        
        
        sJ = imfilter(J,fspecial('disk',50),'replicate');
        tmpHSV = rgb2hsv(sJ);
        tmpU = mean(tmpHSV(:,:,2),2);
        tmpS = std(tmpHSV(:,:,2),1,2);
        tmpSIG = tmpU.*tmpS;
        tmpSIG = bindVec(tmpSIG);
        tmpBV = tmpSIG > graythresh(tmpSIG);
        SUB = sum(tmpBV);
        CB{b}(4) = CB{b}(4) - SUB;
        %{
        imshow(J,[]);
        drawnow
        %}
    end
    
    
    
    
    for e = 1:numel(CB)
        tmpI = imcrop(oI,CB{e});
        
        [containerMask] = getPlantContainerMask(tmpI,cluster_Level0,cluster_Level1);
        cutVec = sum(containerMask,2);
        cidx = find(cutVec > 50);
        if ~isempty(cidx)
            % find the amount to cut from bottom
            toCut = size(tmpI,1) - min(cidx);
            % reshape crop box
            CB{e}(4) = CB{e}(4) - toCut;
            % recut image
            tmpI = imcrop(oI,CB{e});
        end




        [mask] = generatePlantMasks_ver2(tmpI,cluster_Level0,cluster_Level1,cluster_Level2);
        % added at phenome 2018 - plant mask was incomplete - floating plant
        bottomIDX = find(any(mask,2));
        if ~isempty(bottomIDX)
            if bottomIDX(end) ~= size(tmpI,1)
                % find the amount to cut from bottom
                toCut = size(tmpI,1) - max(bottomIDX);
                % reshape crop box
                CB{e}(4) = CB{e}(4) - (toCut+1);
                % recut image
                tmpI = imcrop(oI,CB{e});
                if toCut ~= 0
                    % recut the mask
                    mask(end-toCut:end,:) = [];
                end
            end
        end
        
        

        if ~isempty(cluster_Level3)

            mask = gatherPlantnHoods(tmpI,mask,[0 21],kE,kU,cluster_Level3,[1:3],[]);
            soilMask = mask == 2;
            sidx = find(any(soilMask,2));
            sidx(sidx < (size(soilMask,1) - 50)) = [];
            miniBUMP = 3;
            %{
            [sr sc] = find(soilMask);
            sidx = sr == size(soilMask,1);


            soilMask = bwareaopen(soilMask,50);
            soilMask = imfill(soilMask,'holes');
            soilMask = imdilate(soilMask,strel('disk',11,0));
            ridx = find(soilMask);
            mask = mask ~= 2 & mask ~= 0;
            mask(ridx) = 0;

            cutVec = any(mask,2);
            %}
            %cidx = find(cutVec);
            mask = mask ~= 2 & mask ~= 0;
            mask = imfill(mask,'holes');
            cidx = sidx;
            if ~isempty(cidx)
                toCut = size(tmpI,1) - (min(cidx) - miniBUMP);
                mask = mask(1:(min(cidx) - miniBUMP-1),:);
                CB{e}(4) = CB{e}(4) - (toCut+1);
                tmpI = imcrop(oI,CB{e});
            else
                toCut = 0;
            end
            
            
            % added at phenome 2018 - plant mask was incomplete - floating plant
            bottomIDX = find(any(mask,2));
            if ~isempty(bottomIDX)
                if max(bottomIDX) ~= size(tmpI,1)
                    % find the amount to cut from bottom
                    toCut = size(tmpI,1) - max(bottomIDX);
                    % reshape crop box
                    CB{e}(4) = CB{e}(4) - (toCut+1);
                    % recut image
                    tmpI = imcrop(oI,CB{e});
                    if toCut ~= 0
                        % recut the mask
                        mask((end-toCut):end,:) = [];
                    end
                end
            end
            

            %{
            soilNumbers = 1;
            soilMask = mask == 2 | mask == 4;
            bigSoil = bwareaopen(soilMask,50);
            soilHoles = imfill(bigSoil,'holes');
            soilHoles = soilHoles == 1 & bigSoil == 0;
            mask(find(soilHoles)) = soilNumbers;
            %mask = mask ~= 2 & mask ~=0 & mask ~=1;
            mask = mask ~=4 & mask ~=0 & mask ~=2;
            mask = imfill(mask,'holes');
            mask = bwareaopen(mask,50);
            mask = imclose(mask,strel('disk',2));
            cutVec = any(mask,2);
            cidx = find(cutVec);
            if ~isempty(cidx)
                toCut = size(tmpI,1) - max(cidx);
                mask = mask(1:max(cidx),:);
                CB{e}(4) = CB{e}(4) - toCut;
            else
                toCut = 0;
            end
            %}
            %{
            stmpI = imfilter(tmpI,fspecial('average',[1 30]),'replicate');

            % last stage segment for soil removal
            [mask] = applyClusterStage(stmpI,mask,cluster_Level3);
            soilNumbers = 1;
            soilMask = mask == soilNumbers;
            bigSoil = bwareaopen(soilMask,30);
            soilHoles = imfill(bigSoil,'holes');
            soilHoles = soilHoles == 1 & bigSoil == 0;
            mask(find(soilHoles)) = soilNumbers;
            %mask = mask ~= 2 & mask ~=0 & mask ~=1;
            mask = mask ~=4 & mask ~=0;
            mask = imfill(mask,'holes');
            mask = bwareaopen(mask,50);
            mask = imclose(mask,strel('disk',2));
            cutVec = any(mask,2);
            cidx = find(cutVec);
            if ~isempty(cidx)
                toCut = size(tmpI,1) - max(cidx);
                mask = mask(1:max(cidx),:);
                CB{e}(4) = CB{e}(4) - toCut;
            else
                toCut = 0;
            end
        end
        %}
        %}
        % might not need this last stage it the black cut works well






        %{
        % last stage segment
        [mask] = applyClusterStage(tmpI,mask,cluster_Level3);
        mask = mask ~= 2 & mask ~=0;
        mask = imfill(mask,'holes');
        mask = bwareaopen(mask,50);
        cutVec = any(mask,2);
        cidx = find(cutVec);
        if ~isempty(cidx)
            toCut = size(tmpI,1) - max(cidx);
            mask = mask(1:max(cidx),:);
            CB{e}(4) = CB{e}(4) - toCut;
        else
            toCut = 0;
        end
        %}
        mask = connectPlant(mask);

        boundingBoxes{e} = CB{e};
        tmpI = imcrop(oI,CB{e});
        
        
        if ~isempty(ccm)
            tmpI = rgb2lin(uint8(tmpI*255));
            cTMP = [];
            for k = 1:size(tmpI,3)
                tempy = tmpI(:,:,k);
                cTMP = [cTMP tempy(:)];
            end
            cTMP = [cTMP ones(size(cTMP,1),1)];
            cTMP = double(cTMP)*ccm;
            tmpI = reshape(cTMP,size(tmpI));
            tmpI = tmpI/255;
        end
        subI{e} = tmpI;
        MASK{e} = mask;
        connectedMASK{e} = connectPlant(mask);
        SKELETON{e} = getPlantSkelton(mask);
        
    end
    %{
    R = regionprops(mask,'BoundingBox');
    
    for e = 1:numel(R)
      
    end
    
    %}
    %{
    
     tic;
        fprintf(['starting:: get mask \n']);
        % get the mask
        MASK{e} = getMASK_ver0(I{e});
        % connect the mask
        connectedMASK{e} = connectPlant(MASK{e});
        % get the largest object in the mask
        connectedMASK{e} = bwlarge(connectedMASK{e});
        fprintf(['ending:: get mask::' num2str(toc) '\n']);
        tic
        fprintf(['starting:: get skeleton \n']);
        SKELETON{e} = getPlantSkelton(MASK{e});
        fprintf(['ending:: get skeleton::' num2str(toc) '\n']);
%}
end










