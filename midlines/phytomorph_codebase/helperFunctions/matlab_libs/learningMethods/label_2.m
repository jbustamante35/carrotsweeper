function [l] = label_2(data,para)    
    % histogram split    
    l = kmeans(data,para.value,'emptyaction','singleton');
end