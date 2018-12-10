function [Xtrain,Ytrain,YtrainSequence,YtrainSequence2] = stackSlidingData(imageStack,maskStack,windowSize)
    fprintf(['************************************\n']);
    fprintf(['start stacking data \n']);
    fprintf(['************************************\n']);
    cnt = 1;
    Xtrain = {};
    Ytrain = [];
    YtrainSequence = {};
    for e = 1:size(imageStack,4)
        fprintf(['start stacking set:' num2str(e) ':' num2str(size(imageStack,4)) '\n']);
        tmpY = [];
        for c = 1:size(imageStack,2)


            vec = imageStack(:,c,:,e);
            vec = squeeze(vec);
            if ~isempty(maskStack)
                mvec = maskStack(:,c,1,e);
                mvec = logical(mvec);
                label = any(mvec,1);
            end
            
            
            tmpV = [];
            for k = 1:size(vec,2)
                tmpV = [tmpV ; im2col(vec(:,k),[windowSize 1],'sliding')];
            end
            vec = tmpV;


            Xtrain{cnt} = vec*255^-1;
            if ~isempty(maskStack)
                Ytrain = [Ytrain categorical(double(label))];
                YtrainSequence{cnt} = categorical(double(mvec))';
                tmpY = [tmpY label];
            end
            
            
            cnt = cnt + 1;
        end
        YtrainSequence2{e} = tmpY;
        fprintf(['end stacking set:' num2str(e) ':' num2str(size(imageStack,4)) '\n']);
    end
    fprintf(['************************************\n']);
    fprintf(['stop stacking data \n']);
    fprintf(['************************************\n']);
end