function [I] = normalize(I)    
    if ~strcmp(class(I),'double')
        I = double(I);
    end
    I = I - min(I(:));
    I = I / max(I(:));
end