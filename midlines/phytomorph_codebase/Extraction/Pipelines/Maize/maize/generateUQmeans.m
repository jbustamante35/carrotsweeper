function [udata UQ] =  generateUQmeans(data,labels)
    UQ = unique(labels);
    for u = 1:numel(UQ)
        udata(u,:) = mean(data(strcmp(labels,UQ{u}),:),1);
    end    
end