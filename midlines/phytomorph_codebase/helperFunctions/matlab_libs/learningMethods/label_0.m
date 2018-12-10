function [l] = label_0(data,para)    
    % histogram split
    level = graythresh(data);
    l = data > level;
end