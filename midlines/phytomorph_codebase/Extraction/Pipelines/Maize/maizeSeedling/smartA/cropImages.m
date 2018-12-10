function [returnI,boundingBoxes] = cropImages(I,reSize,networkObject1,networkObject2,networkObject3)
    fprintf(['..................................\n']);
    fprintf(['..................................\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % needed constants
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SNIP = 10;
    FULL = 200;
    PAD = 5;
    toCLOSE = 700;
    WID_INCREASE = .35;
    BOT_THRESH = 300;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read image if file name passed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ischar(I)
        tic;
        fprintf(['start::reading file \n']);
        I = double(imread(I));
        [I angle] = rectifyImage(I);
        fprintf(['end::reading file::' num2str(toc) '\n']);
    end
    %I(:,1:70,:) = [];
    %I(:,(end-69):end,:) = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read image if file name passed
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % do vertical strip sizing
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %oSZ = size(I);
    %I = imresize(I,[3280 size(I,2)]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % obtain vertical strips
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    fprintf(['start::finding vertical strips \n']);
    % resize permute and reshape
    fprintf(['start::resize-permute-reshape \n']);
    rI = imresize(I,reSize);
    rI = permute(rI,[1 3 2]);
    rI = reshape(rI,[size(rI,1)*size(rI,2) size(rI,3)]);
    fprintf(['end::resize-permute-reshape \n']);
    % apply neural network
    fprintf(['start::apply NN \n']);
    sig1 = applyScannerNetwork(networkObject1,rI);
    sig1 = imclose(imresize(sig1 > .5,[1 size(I,2)]),ones([1 50]));
    C = repmat(sig1,[size(I,1) 1]);
    fprintf(['start::end NN \n']);
    % crop vertical strips
    fprintf(['start::cropping vertical strips \n']);
    [sI,~,realS,R] = verticalCrop_ver0(I,C,1200);
    fprintf(['end::cropping vertical strips \n']);
    fprintf(['end::finding vertical strips::' num2str(toc) '\n']);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % obtain vertical strips
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TOP_THRESH = 950;
    fprintf(['..................................\n']);
    fprintf(['..................................\n']);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % find horizontal cropping for each binary object
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for e = 1:size(sI,4)
        fprintf(['..................................\n']);
        fprintf(['..................................\n']);
        tic;
        fprintf(['start::final cropping \n']);
        % make mask
        fprintf(['start::making mask \n']);
        Z = zeros([size(I,1) size(I,2)]);
        fprintf(['end::making mask \n']);
        % resize permute and reshape
        fprintf(['start::resize-permute-reshape \n']);
        tmp = sI(:,:,:,e);
        tmp = sort(tmp,2);
        tmp = permute(tmp,[2 1 3]);
        tmp = permute(tmp,[1 3 2]);
        tmp = reshape(tmp,[size(tmp,1)*size(tmp,2) size(tmp,3)]);
        fprintf(['end::resize-permute-reshape \n']);
        % apply neural network
        fprintf(['start::apply NN \n']);
        sig2 = applyScannerNetwork(networkObject2,tmp);
        sig2(1:TOP_THRESH) = 0;
        sig2((end-BOT_THRESH):end) = 0;
        sig2B = sig2 > .2;
        sig2B = bwareaopen(sig2B,70);
        sig2 = imclose(sig2B,ones([1 toCLOSE]))';
        idx = find(sig2);
        sig2(idx(end-SNIP:end)) = 0;
        fprintf(['end::apply NN \n']);
        fprintf(['start::preparing mask \n']);
        % get bounding box
        BOX = round(R(e).BoundingBox);
        % fill in mask
        Z(BOX(2):(BOX(2)+BOX(4)-1),BOX(1):(BOX(1)+BOX(3)-1)) = repmat(sig2,[1 (BOX(1)+BOX(3)-BOX(1))]);
        % get largest object
        Z = bwlarge(Z);
        % find object
        tmpR = regionprops(logical(Z));
        fprintf(['end::preparing mask \n']);
        % crop the vertical and horizontal
        fprintf(['start::cropping vertical+horizontal \n']);
        subI = imcrop(I,tmpR(1).BoundingBox);
        fprintf(['end::cropping vertical+horizontal \n']);
        % resize permute and reshape
        fprintf(['start::resize-permute-reshape \n']);
        bot = imresize(subI(end-FULL:end,:,:),[FULL+1 1200]);
        bot = sort(bot,2);
        bot = permute(bot,[2 1 3]);
        bot = permute(bot,[1 3 2]);
        bot = reshape(bot,[size(bot,1)*size(bot,2) size(bot,3)]);
        fprintf(['end::resize-permute-reshape \n']);
        % apply neural network
        fprintf(['start::apply NN \n']);
        sig3 = applyScannerNetwork(networkObject3,bot);
        sig3 = sum(sig3 > .7) + PAD;
        fprintf(['end::apply NN \n']);
        % prepare final bounding box
        fprintf(['start::preparing final bounding box \n']);
        toAdd = round(tmpR(1).BoundingBox(3)*WID_INCREASE/2);
        % make the box wider 
        tmpR(1).BoundingBox(1) = max(tmpR(1).BoundingBox(1)-toAdd,1);
        tmpR(1).BoundingBox(3) = tmpR(1).BoundingBox(3) + 2*toAdd;
        % make it shorter via move
        tmpR(1).BoundingBox(4) = tmpR(1).BoundingBox(4) - sig3 - 70;
        tmpR(1).BoundingBox(2) = tmpR(1).BoundingBox(2)+70;
        % make it a bit taller via stretch
        ST = 40;
        tmpR(1).BoundingBox(4) = tmpR(1).BoundingBox(4) + ST;
        tmpR(1).BoundingBox(2) = tmpR(1).BoundingBox(2) - ST;
        fprintf(['start::preparing final bounding box \n']);
        % crop the image at the final bounding box
        fprintf(['start::final crop for sub image \n']);
        subI = imcrop(I,tmpR(1).BoundingBox);
        fprintf(['end::final crop for sub image \n']);
        % store the image for return
        returnI{e} = subI;
        boundingBoxes{e} = tmpR(1).BoundingBox;
        fprintf(['start::final cropping::' num2str(toc) '\n']);
        fprintf(['..................................\n']);
        fprintf(['..................................\n']);
    end
end