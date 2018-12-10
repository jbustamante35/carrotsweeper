function [Xtrain,Ytrain,numClasses] = stackColumnData(imageStack,maskStack,toSort)
    fprintf(['************************************\n']);
    fprintf(['start stacking data \n']);
    fprintf(['************************************\n']);
    cnt = 1;
    Xtrain = {};
    Ytrain = {};
    YtrainSequence = {};
    numClasses = 0;
    for e = 1:size(imageStack,4)
        fprintf(['start stacking set:' num2str(e) ':' num2str(size(imageStack,4)) '\n']);
        tmpY = [];
        tmpX = [];
        for c = 1:size(imageStack,2)

            vec = imageStack(:,c,:,e);
            vec = squeeze(vec);
            if toSort
                vec = sort(vec,1);
            end
            
            if ~isempty(maskStack)
                mvec = maskStack(:,c,1,e);
                mvec = logical(mvec);
                label = any(mvec,1);
            end
            
            
            
           tmpX =[tmpX vec(:)*255^-1];
            
            
            if ~isempty(maskStack)
                tmpY = [tmpY label];
            end
            
            
            cnt = cnt + 1;
        end
        Xtrain{e} = tmpX;
        
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
        Ytrain{e} = categorical(newLabel);
        
        fprintf(['end stacking set:' num2str(e) ':' num2str(size(imageStack,4)) '\n']);
    end
    fprintf(['************************************\n']);
    fprintf(['stop stacking data \n']);
    fprintf(['************************************\n']);
end