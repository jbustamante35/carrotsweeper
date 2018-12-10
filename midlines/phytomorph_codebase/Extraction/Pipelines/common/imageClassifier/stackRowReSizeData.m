function [Xtrain,Ytrain,numClasses] = stackRowReSizeData(imageStack,maskStack,Resize)
    fprintf(['************************************\n']);
    fprintf(['start stacking data \n']);
    fprintf(['************************************\n']);
    cnt = 1;
    Xtrain = {};
    Ytrain = {};
    YtrainSequence = {};
    for e = 1:size(imageStack,4)
        fprintf(['start stacking set:' num2str(e) ':' num2str(size(imageStack,4)) '\n']);
        tmpY = [];
        tmpX = [];
        
        R = regionprops(maskStack(:,:,1,e)>=1,'BoundingBox');
        
        numClasses = 0;
        
        for r = 1:numel(R)
            recROW = 1:size(imageStack,1);
            recCOL = round(R(r).BoundingBox(1):(R(r).BoundingBox(1) + R(r).BoundingBox(3)));
            recIMG = imageStack(recROW,recCOL,:,e);
            
            if ~isempty(maskStack)
                recMSK = maskStack(recROW,recCOL,1,e);
                tmpY = imresize(recMSK,[size(recMSK,1) Resize]);
                tmpY = (any(tmpY>=1,2));
                tmpY = permute(tmpY,[2 1]);
                
                
                
                
                % relabel sequence
                label = tmpY == 1;
                mR0 = regionprops(label==0,'PixelIdxList');
                mR1 = regionprops(label==1,'PixelIdxList');
                newLabel = zeros(size(label));
                for r = 1:numel(mR0)
                    value = 2*(r-1);
                    newLabel(mR0(r).PixelIdxList) = value;
                end
                for r = 1:numel(mR1)
                    value = 2*(r-1)+1;
                    newLabel(mR1(r).PixelIdxList) = value;
                end
                
                numClasses = max(numClasses,numel(unique(newLabel)));
                
                Ytrain{cnt} = categorical(newLabel);
            end
            
            tmpX = imresize(recIMG,[size(recIMG,1) Resize]);
            
            
            tmpX = permute(tmpX,[2 1 3]);
            tmpX = permute(tmpX,[1 3 2]);
            sz = size(tmpX);
            
            
            Xtrain{cnt} = reshape(tmpX,[prod(sz(1:2)) sz(3)]);
            
            
            
            cnt = cnt + 1;
        end
        
        
        fprintf(['end stacking set:' num2str(e) ':' num2str(size(imageStack,4)) '\n']);
    end
    fprintf(['************************************\n']);
    fprintf(['stop stacking data \n']);
    fprintf(['************************************\n']);
end