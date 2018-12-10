function [box] = mergeCropBoxes(box)
    boxStack = [];
    for e = 1:numel(box)
        tmp = [];
        for b = 1:numel(box{e})
            tmp = [tmp;box{e}{b}];
        end
        boxStack = cat(3,boxStack,tmp);
    end
    boxStack = mean(boxStack,3);
    
    meanDIM = mean(boxStack(:,3:4),1);
    boxStack(:,3) = meanDIM(1);
    boxStack(:,4) = meanDIM(2);
    
    boxStack(:,1) = round(boxStack(:,1));
    boxStack(:,2) = round(boxStack(:,2));
    
    box = [];
    for e = 1:size(boxStack,1)
        box{e} = boxStack(e,:);
    end
end