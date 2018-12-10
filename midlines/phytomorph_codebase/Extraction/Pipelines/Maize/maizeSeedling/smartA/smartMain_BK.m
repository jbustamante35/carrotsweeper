function [measureBlocks] = smartMain(fileName,reSize,networkObject1,networkObject2,networkObject3,oPath,rPath,toM)

    fprintf(['************************************ STARTING STAGE ONE ************************************ \n']);
    fprintf(['************************************    image crop      ************************************ \n']);
    % crop the plants with the NN
    I = double(imread(fileName));
    [returnI,boundingBoxes] = cropImages(I,reSize,networkObject1,networkObject2,networkObject3);
    fprintf(['************************************  ENDING STAGE ONE  ************************************ \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['************************************ STARTING STAGE TWO ************************************ \n']);
    fprintf(['************************************  mask and skeleton ************************************ \n']);
    % get the mask(s) and skeleton(s)
    [connectedMASK,MASK,SKELETON] = getMASKandSKELETON(returnI);
    fprintf(['************************************  ENDING STAGE TWO  ************************************ \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['************************************STARTING STAGE THREE************************************ \n']);
    fprintf(['************************************ measure phenotypes ************************************ \n']);
    [phenoTypes measureBlocks] = measurePhenotypes(MASK,SKELETON,toM);
    fprintf(['************************************ ENDING STAGE THREE ************************************ \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['******************************************************************************************** \n']);
    fprintf(['************************************STARTING STAGE FOUR ************************************ \n']);
    fprintf(['************************************  display and save  ************************************ \n']);
    displayResults(fileName,I,returnI,boundingBoxes,MASK,SKELETON,phenoTypes,oPath,rPath);
    fprintf(['************************************  ENDING STAGE FOUR ************************************ \n']);
    
    %{
    extract = 0;
    E = [];
    I = double(imread(fileName));
    rI = imresize(I,reSize);
    rI = permute(rI,[1 3 2]);
    rI = reshape(rI,[size(rI,1)*size(rI,2) size(rI,3)]);
    
    % apply network
    sig1 = applyScannerNetwork(networkObject1,rI);
    sig1 = imclose(imresize(sig1 > .5,[1 size(I,2)]),ones([1 50]));
    C = repmat(sig1,[size(I,1) 1]);
    
    [sI,~,realS,R] = verticalCrop_ver0(I,C,1200);
    
    SNIP = 10;
    FULL = 200;
    PAD = 3;
    WID_INCREASE = .35;
    for e = 1:size(sI,4)
        Z = zeros([size(I,1) size(I,2)]);
        
        tmp = sI(:,:,:,e);
        tmp = permute(tmp,[2 1 3]);
        tmp = permute(tmp,[1 3 2]);
        tmp = reshape(tmp,[size(tmp,1)*size(tmp,2) size(tmp,3)]);
        sig2 = applyScannerNetwork(networkObject2,tmp);
        
        sig2 = imclose(sig2 > .7,ones([1 200]))';
        idx = find(sig2);
        sig2(idx(end-SNIP:end)) = 0;
        
        BOX = round(R(e).BoundingBox);
        
        
        Z(BOX(2):(BOX(2)+BOX(4)-1),BOX(1):(BOX(1)+BOX(3)-1)) = repmat(sig2,[1 (BOX(1)+BOX(3)-BOX(1))]);
        Z = bwlarge(Z);
        
        
        tmpR = regionprops(logical(Z));
        subI = imcrop(I,tmpR(1).BoundingBox);
        
        bot = imresize(subI(end-FULL:end,:,:),[FULL+1 1200]);
        bot = permute(bot,[2 1 3]);
        bot = permute(bot,[1 3 2]);
        bot = reshape(bot,[size(bot,1)*size(bot,2) size(bot,3)]);
        sig3 = applyScannerNetwork(networkObject3,bot);
        sig3 = sum(sig3 > .7) + PAD;
        
        toAdd = round(tmpR(1).BoundingBox(1)*WID_INCREASE/2);
        % make the box wider 
        tmpR(1).BoundingBox(1) = max(tmpR(1).BoundingBox(1)-toAdd,1);
        tmpR(1).BoundingBox(3) = tmpR(1).BoundingBox(3) + toAdd;
        
        % make it shorter
        tmpR(1).BoundingBox(4) = tmpR(1).BoundingBox(4) - sig3 - 50;
        tmpR(1).BoundingBox(2) = tmpR(1).BoundingBox(2) + 50;
        subI = imcrop(I,tmpR(1).BoundingBox);
        
        
        
        
        
        fI{e} = subI;
        
        %}
        %{
        STEM_SNIP = 20;
        outLN = 4;
        MASK = getMASK_ver0(subI);
        MASK = connectPlant(MASK);
        MASK = bwlarge(MASK);
        [skeleton] = getPlantSkelton(MASK);
        [basePoint] = getStemBasePoint(MASK,SNIP,skeleton);
        [diameter] = measureStemDiameter(MASK,STEM_SNIP,outLN,basePoint);
        
        
        ************************************************************************%%
        % trace the skeleton
        ************************************************************************%%
        [path,pathcost,midx,DP] = traceSeedlingSkeleton(skeleton,basePoint);
        ************************************************************************%%
        % trace the skeleton
        ************************************************************************%%
        
        
        
        plot(CenterOfMass(2),CenterOfMass(1),'c*');
        plot(linspace((CenterOfMass(2)-WIDTH(e)/2),(CenterOfMass(2)+WIDTH(e)/2),2),[CenterOfMass(1) CenterOfMass(1)],'c');
        plot([CenterOfMass(2) CenterOfMass(2)],linspace((CenterOfMass(1)-plantHEIGHT(e)/2),(CenterOfMass(1)+plantHEIGHT(e)/2),2),'c');
        TH = linspace(-pi,pi,200);
        X = StdDis(2)*cos(TH) + CenterOfMass(2);
        Y = StdDis(1)*sin(TH) + CenterOfMass(1);
        
        if extract
            E(:,:,:,e) = imresize(subI(end-FULL:end,:,:),[FULL+1 1200]);
        end
        
    end
    %}
    
end