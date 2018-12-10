function [samp samp_Mask topIdx botIdx] = sampleBoundary(image,dB,NOR,scale,tip,alongValue,disp)
    samp = [];
    image = double(imread(image));
    image = image(:,:,1:3);
    image = permute(image,[2 1 3]);
    G = rgb2gray(image/255);
    dG = stdfilt(G,ones(11));
    
    dL = diff(dB{1},1,1);
    dL = sum(dL.*dL,2).^.5;
    L = cumsum([0;dL],1);
    tL = L(tip);
    topIdx = find((L - tL) + alongValue > 0);
    botIdx = find((L - tL) - alongValue < 0);
    
    tmpSampBot = [];
    tmpSampTop = [];
    tmpSampTop_Mask = [];
    tmpSampBot_Mask = [];
    for s = scale
        ED = [];
        for p = 1:size(dB{1},1)
            str = dB{1}(p,:);
            ed = str +s*NOR(p,:);            
            ED = [ED;ed];
        end
        
        tmpSampTop = cat(3,tmpSampTop,squeeze(ba_interp2(image,ED(topIdx(1):tip,2),ED(topIdx(1):tip,1))));
        tmpSampTop_Mask = cat(2,tmpSampTop_Mask,squeeze(ba_interp2(dG,ED(topIdx(1):tip,2),ED(topIdx(1):tip,1))));
        
        tmpSampBot = cat(3,tmpSampBot,squeeze(ba_interp2(image,ED(tip:botIdx(end),2),ED(tip:botIdx(end),1),1)));
        tmpSampBot_Mask = cat(2,tmpSampBot_Mask,squeeze(ba_interp2(dG,ED(tip:botIdx(end),2),ED(tip:botIdx(end),1),1)));
    end
    
    % create mask for the top and average along the "radial" direction
    %level = graythresh(tmpSampTop_Mask);
    %tmpSampTop_Mask = tmpSampTop_Mask > level;
    for e = 1:size(tmpSampTop_Mask,1)
        level = graythresh(tmpSampTop_Mask(e,:));
        tmpSampTop_Mask(e,:) = tmpSampTop_Mask(e,:) > level;
    end
    
    zerosMaskTop = zeros(size(tmpSampTop_Mask));
    for e = 1:size(tmpSampTop_Mask,1)
        tmp = tmpSampTop_Mask(e,:);
        fidx = find(tmp);
        if ~isempty(fidx)
            idx = fidx(end);
        end
        stop = min(idx+7,size(zerosMaskTop,2));
        zerosMaskTop(e,idx:stop) = 1;
    end
    tmpSampTop = permute(tmpSampTop,[1 3 2]);
    tmpSampTop = tmpSampTop.*repmat(zerosMaskTop,[1 1 3]);
    tmpSampTop = squeeze(sum(tmpSampTop,2)).*repmat(sum(zerosMaskTop,2),[1 3]).^-1;
    tmpSampTop = interp1(1:size(tmpSampTop,1),tmpSampTop,linspace(1,size(tmpSampTop,1),alongValue));
    
    % create mask for the bottom and average along the "radial" direction
    %level = graythresh(tmpSampBot_Mask);
    %tmpSampBot_Mask = tmpSampBot_Mask > level;
    for e = 1:size(tmpSampBot_Mask,1)
        level = graythresh(tmpSampBot_Mask(e,:));
        tmpSampBot_Mask(e,:) = tmpSampBot_Mask(e,:) > level;
    end
    
    zerosMaskBot = zeros(size(tmpSampBot_Mask));
    for e = 1:size(tmpSampBot_Mask,1)
        tmp = tmpSampBot_Mask(e,:);
        fidx = find(tmp);
        if ~isempty(fidx)
            idx = fidx(end);
        end
        stop = min(idx+7,size(zerosMaskBot,2));
        zerosMaskBot(e,idx:stop) = 1;
    end
    tmpSampBot = permute(tmpSampBot,[1 3 2]);
    tmpSampBot = tmpSampBot.*repmat(zerosMaskBot,[1 1 3]);
    tmpSampBot = squeeze(sum(tmpSampBot,2)).*repmat(sum(zerosMaskBot,2),[1 3]).^-1;
    tmpSampBot = interp1(1:size(tmpSampBot,1),tmpSampBot,linspace(1,size(tmpSampBot,1),alongValue));

    
    %{
    tmpSampBot = permute(tmpSampBot,[1 3 2]);
    tmpSampBot = squeeze(mean(tmpSampBot,2));
    tmpSampBot = interp1(1:size(tmpSampBot,1),tmpSampBot,linspace(1,size(tmpSampBot,1),alongValue));
    %}
    
    
    samp = cat(1,tmpSampTop,tmpSampBot);
    samp_Mask = cat(1,tmpSampTop_Mask,tmpSampBot_Mask);
    
    
    
    
    
    if disp
        imshow(image)
        hold on
        plot(dB{1}(topIdx(1):tip,2),dB{1}(topIdx(1):tip,1),'r');
        plot(dB{1}(tip:botIdx(end),2),dB{1}(tip:botIdx(end),1),'b');
        drawnow
    end
end